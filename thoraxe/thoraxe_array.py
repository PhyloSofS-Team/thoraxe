"""
thoraxe_array: run ThorAxe on every gene folder of a directory, with a SLURM job array.

Usage, on a login node, with the Python where thoraxe is installed::

    python thoraxe/thoraxe_array.py [options] IN_DIR OUT_DIR -- [thoraxe options]
    python thoraxe/thoraxe_array.py --resume OUT_DIR    # submit again; finished genes are skipped
    python thoraxe/thoraxe_array.py --status OUT_DIR    # progress and failed genes

IN_DIR holds one folder per gene with an Ensembl/ subfolder, as written by transcript_query
(or transcript_query_bulk.sbatch), with or without add_transcripts. ThorAxe runs on every
gene folder whose Ensembl/exonstable.tsv exists (genes without orthologs have none) and
writes OUT_DIR/<gene>/thoraxe/, as ``thoraxe -i IN_DIR/<gene> -o OUT_DIR/<gene>`` followed
by the thoraxe options given after ``--``.

Before submitting, OUT_DIR/run/ receives what the run depends on, so that it can be
reproduced and is not affected by later changes to the ThorAxe clone:

- ``code/``: the thoraxe package of the installed git commit (tracked files only;
  uncommitted changes are refused);
- ``bin/ProGraphMSA``: a copy of the aligner;
- ``species_list.txt``: a copy of the -l file, if one is given (the option is rewritten
  to point at the copy);
- ``genes.txt``: the genes to run, largest exon table first, so that the slow genes start
  first and land in different array tasks; ``genes_without_input.txt``: the others;
- ``thoraxe_args.txt``, ``input_dir``, ``python``, ``seed``: read by the worker,
  ``thoraxe_array.sbatch``, copied there and submitted;
- ``pip_freeze.txt``: the Python packages.

OUT_DIR/thoraxe_info.json describes the run: input, ThorAxe commit and every parameter
(defaults included), species, aligner and Python versions, checksums, and each submission.
``--provenance FILE`` adds a JSON description of the input (how IN_DIR was built).

Each gene is logged in OUT_DIR/thoraxe_logs/<gene>.log, which ends with its wall time and
peak memory and ``exit=<code>``; SLURM logs go to OUT_DIR/slurm_logs/. The sbatch options
default to the worker's (``--array=0-999``, 1 CPU, 8 GB, 1 day per task); ``--sbatch``
adds or overrides some, e.g. ``--sbatch "--qos=std_long --time=3-00:00:00"``.
"""

import argparse
import datetime
import getpass
import hashlib
import json
import os
import shlex
import shutil
import socket
import subprocess
import sys
import tarfile
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WORKER = os.path.join(HERE, 'thoraxe_array.sbatch')
INFO = 'thoraxe_info.json'

# Import the thoraxe of the snapshot given as first argument (sys.path[0]: an editable
# install would win over PYTHONPATH), check it is the one imported, then print the version
# and every parsed option (defaults included). A bad option exits with thoraxe's usage.
_PARSE = r'''
import json, sys
code = sys.argv.pop(1)
sys.path.insert(0, code)
import thoraxe
if not thoraxe.__file__.startswith(code):
    sys.exit("thoraxe imported from %s, not from %s" % (thoraxe.__file__, code))
from thoraxe.thoraxe import parse_command_line
from thoraxe.version import __version__
args = vars(parse_command_line().parse_args(sys.argv[1:]))
args.pop("inputdir")
args.pop("outputdir")
print(json.dumps({"version": __version__, "parameters": args}))
'''

_VERSIONS = r'''
import json, platform
out = {"python": platform.python_version()}
for name in ("Bio", "numpy", "pandas", "scipy", "networkx", "plotly"):
    try:
        out[name] = __import__(name).__version__
    except Exception:  # pylint: disable=broad-except
        out[name] = None
print(json.dumps(out))
'''


def _now():
    return datetime.datetime.now().astimezone().isoformat(timespec='seconds')


def _run(cmd, env=None):
    "Return the stripped standard output of cmd; raise with its error output if it fails."
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, env=env, check=False)
    if proc.returncode != 0:
        raise SystemExit(f'failed: {" ".join(cmd)}\n{proc.stderr.strip()}')
    return proc.stdout.strip()


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _clean_env():
    env = dict(os.environ)
    env.pop('PYTHONPATH', None)
    env['OPENBLAS_NUM_THREADS'] = '1'
    return env


def _git(repo, *args):
    return _run(['git', '-C', repo] + list(args))


def thoraxe_repo(python):
    "Return the git clone that python imports thoraxe from, with its commit data."
    package = _run([python, '-c', 'import thoraxe; print(thoraxe.__file__)'], env=_clean_env())
    repo = os.path.dirname(os.path.dirname(os.path.abspath(package)))
    if not os.path.isdir(os.path.join(repo, '.git')):
        raise SystemExit(f'thoraxe is imported from {package}, which is not in a git clone: '
                         'the run must record a commit')
    dirty = _git(repo, 'status', '--porcelain', '--untracked-files=no', '--', 'thoraxe')
    if dirty:
        raise SystemExit(f'{repo} has uncommitted changes in thoraxe/:\n{dirty}\n'
                         'Commit or stash them: the run records the commit it uses.')
    commit = _git(repo, 'rev-parse', 'HEAD')
    try:
        remote = _git(repo, 'remote', 'get-url', 'origin')
    except SystemExit:
        remote = None
    return repo, {
        'repository': repo,
        'commit': commit,
        'describe': _git(repo, 'describe', '--tags', '--always'),
        'branch': _git(repo, 'rev-parse', '--abbrev-ref', 'HEAD'),
        'commit_date': _git(repo, 'show', '-s', '--format=%cI', commit),
        'remote': remote,
        'on_remote_branch': bool(_git(repo, 'branch', '-r', '--contains', commit)),
    }


def snapshot(repo, dest):
    "Extract the thoraxe package of HEAD into dest."
    os.makedirs(dest)
    with subprocess.Popen(['git', '-C', repo, 'archive', '--format=tar', 'HEAD', 'thoraxe'],
                          stdout=subprocess.PIPE) as proc:
        with tarfile.open(fileobj=proc.stdout, mode='r|') as archive:
            archive.extractall(dest)
    if proc.returncode != 0:
        raise SystemExit(f'git archive failed in {repo}')


def _option_value(args, short, long):
    "Return (index, value) of an option written '-x V', '--long V' or '--long=V'."
    for i, arg in enumerate(args):
        if arg in (short, long) and i + 1 < len(args):
            return i, args[i + 1]
        if arg.startswith(long + '='):
            return i, arg.split('=', 1)[1]
    return None, None


def _set_option(args, index, value):
    "Replace the value of the option found at index by _option_value."
    if args[index].startswith('--') and '=' in args[index]:
        args[index] = args[index].split('=', 1)[0] + '=' + value
    else:
        args[index + 1] = value


def list_genes(in_dir, only=None):
    """
    Return (genes to run, largest exon table first; gene folders without exon table).

    A gene folder is a subfolder of in_dir with an Ensembl/ subfolder.
    """
    with_input, without = [], []
    names = sorted(only) if only is not None else sorted(os.listdir(in_dir))
    for name in names:
        ensembl = os.path.join(in_dir, name, 'Ensembl')
        if not os.path.isdir(ensembl):
            if only is not None:
                raise SystemExit(f'{name}: no {ensembl}')
            continue
        table = os.path.join(ensembl, 'exonstable.tsv')
        size = os.path.getsize(table) if os.path.isfile(table) else 0
        if size:
            with_input.append((size, name))
        else:
            without.append(name)
    with_input.sort(key=lambda item: (-item[0], item[1]))
    return [name for _, name in with_input], without


def _species_seen(in_dir, genes, count=20):
    "Species of the tsl.csv files of the first genes."
    seen = set()
    for gene in genes[:count]:
        with open(os.path.join(in_dir, gene, 'Ensembl', 'tsl.csv'), encoding='utf-8') as stream:
            next(stream, None)
            seen.update(line.split(',', 1)[0] for line in stream if line.strip())
    return seen


def _write(path, text):
    with open(path, 'w', encoding='utf-8') as stream:
        stream.write(text)


def _write_info(out_dir, info):
    path = os.path.join(out_dir, INFO)
    _write(path + '.tmp', json.dumps(info, indent=1) + '\n')
    os.replace(path + '.tmp', path)


def prepare(opts, in_dir, out_dir, thoraxe_args, run_dir):  # pylint: disable=too-many-locals,too-many-statements
    "Fill run_dir with the run inputs and return the content of thoraxe_info.json."
    reserved = {'-i', '--inputdir', '-o', '--outputdir'}
    used = [arg for arg in thoraxe_args if arg.split('=', 1)[0] in reserved]
    if used:
        raise SystemExit(f'{" ".join(used)}: input and output folders are set per gene')
    python = os.path.abspath(opts.python)
    repo, git = thoraxe_repo(python)
    code_dir = os.path.join(run_dir, 'code')
    snapshot(repo, code_dir)
    args = list(thoraxe_args)

    # aligner: thoraxe's -a, or ProGraphMSA next to python or on PATH
    index, aligner = _option_value(args, '-a', '--aligner')
    search = os.pathsep.join([os.path.dirname(python), os.environ.get('PATH', '')])
    found = shutil.which(aligner or 'ProGraphMSA', path=search)
    if not found:
        raise SystemExit(f'aligner {aligner or "ProGraphMSA"} not found')
    found = os.path.realpath(found)
    os.makedirs(os.path.join(run_dir, 'bin'))
    aligner_copy = os.path.join(run_dir, 'bin', 'ProGraphMSA')
    shutil.copy2(found, aligner_copy)
    if index is not None:
        _set_option(args, index, aligner_copy)
    aligner_version = subprocess.run([aligner_copy, '--version'], stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, universal_newlines=True,
                                     check=False).stdout.strip()
    aligner_version = aligner_version.replace(aligner_copy, '').strip()

    # genes
    only = None
    if opts.genes:
        with open(opts.genes, encoding='utf-8') as stream:
            only = {line.strip() for line in stream if line.strip()}
    genes, without = list_genes(in_dir, only)
    if not genes:
        raise SystemExit(f'no gene folder with Ensembl/exonstable.tsv in {in_dir}')
    _write(os.path.join(run_dir, 'genes.txt'), ''.join(g + '\n' for g in genes))
    _write(os.path.join(run_dir, 'genes_without_input.txt'), ''.join(g + '\n' for g in without))

    # species list: copied, and checked against the first genes
    species_info = None
    index, species = _option_value(args, '-l', '--specieslist')
    if species:
        if os.path.isfile(species):
            with open(species, encoding='utf-8') as stream:
                names = [line.strip() for line in stream if line.strip()]
            copy = os.path.join(run_dir, 'species_list.txt')
            _write(copy, ''.join(name + '\n' for name in names))
            _set_option(args, index, copy)
            species_info = {'file': 'run/species_list.txt', 'source': os.path.abspath(species),
                            'species': names}
        else:
            names = species.split(',')
            species_info = {'species': names}
        unseen = sorted(set(names) - _species_seen(in_dir, genes))
        if unseen:
            print(f'warning: not in the tsl.csv of the 20 largest genes: {", ".join(unseen)}')

    parsed = json.loads(_run([python, '-c', _PARSE, code_dir + os.sep] + args,
                             env=_clean_env()))
    versions = json.loads(_run([python, '-c', _VERSIONS], env=_clean_env()))
    try:
        freeze = _run([python, '-m', 'pip', 'freeze'], env=_clean_env())
    except SystemExit:
        freeze = ''
    _write(os.path.join(run_dir, 'pip_freeze.txt'), freeze + '\n')

    _write(os.path.join(run_dir, 'thoraxe_args.txt'), ''.join(a + '\n' for a in args))
    _write(os.path.join(run_dir, 'input_dir'), in_dir + '\n')
    _write(os.path.join(run_dir, 'python'), python + '\n')
    _write(os.path.join(run_dir, 'seed'), f'{opts.seed}\n')
    shutil.copy2(WORKER, os.path.join(run_dir, 'thoraxe_array.sbatch'))

    version_csv = os.path.join(in_dir, genes[0], 'Ensembl', 'ensembl_version.csv')
    provenance = None
    if opts.provenance:
        with open(opts.provenance, encoding='utf-8') as stream:
            provenance = json.load(stream)
    return {
        'created': _now(),
        'user': getpass.getuser(),
        'host': socket.gethostname(),
        'input': {
            'dir': in_dir,
            'gene_folders': len(genes) + len(without),
            'genes_run': len(genes),
            'genes_without_exon_table': len(without),
            'genes_file': 'run/genes.txt',
            'genes_file_sha256': _sha256(os.path.join(run_dir, 'genes.txt')),
            'ensembl_version_csv': (open(version_csv, encoding='utf-8').read().strip()  # pylint: disable=consider-using-with
                                    if os.path.isfile(version_csv) else None),
            'provenance': provenance,
        },
        'output': {
            'dir': out_dir,
            'per_gene': '<gene>/thoraxe/',
            'logs': 'thoraxe_logs/<gene>.log (last lines: resources, exit code)',
            'slurm_logs': 'slurm_logs/',
        },
        'thoraxe': {
            'version': parsed['version'],
            'git': git,
            'code_snapshot': 'run/code',
            'command': 'thoraxe -i IN_DIR/<gene> -o OUT_DIR/<gene> ' + ' '.join(
                shlex.quote(a) for a in thoraxe_args),
            'arguments_as_run': args,
            'parameters': parsed['parameters'],
            'environment': {'PYTHONHASHSEED': str(opts.seed), 'OPENBLAS_NUM_THREADS': '1',
                            'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'},
        },
        'species_list': species_info,
        'aligner': {
            'source': found,
            'copy': 'run/bin/ProGraphMSA',
            'sha256': _sha256(aligner_copy),
            'version': aligner_version,
        },
        'python': {
            'executable': python,
            'versions': versions,
            'pip_freeze': 'run/pip_freeze.txt',
        },
        'launcher': {
            'script': os.path.abspath(__file__),
            'sha256': _sha256(os.path.abspath(__file__)),
            'worker': 'run/thoraxe_array.sbatch',
            'worker_sha256': _sha256(WORKER),
        },
        'submissions': [],
    }


def submit(out_dir, sbatch_options, dry_run=False):
    "Submit the worker of out_dir/run; return the job id (None for a dry run)."
    cmd = ['sbatch', '--parsable', f'--output={out_dir}/slurm_logs/%A_%a.log'] + \
        shlex.split(sbatch_options or '') + \
        [os.path.join(out_dir, 'run', 'thoraxe_array.sbatch'), out_dir]
    print('$ ' + ' '.join(shlex.quote(c) for c in cmd))
    if dry_run:
        return None
    os.makedirs(os.path.join(out_dir, 'slurm_logs'), exist_ok=True)
    job = _run(cmd).split(';')[0]
    info_path = os.path.join(out_dir, INFO)
    with open(info_path, encoding='utf-8') as stream:
        info = json.load(stream)
    info['submissions'].append({'job_id': job, 'date': _now(), 'command': cmd})
    _write_info(out_dir, info)
    print(f'submitted job {job}')
    return job


def status(out_dir):
    "Print the progress of out_dir."
    run_dir = os.path.join(out_dir, 'run')
    with open(os.path.join(run_dir, 'genes.txt'), encoding='utf-8') as stream:
        genes = [line.strip() for line in stream if line.strip()]
    codes, failed, wall = {}, [], 0.0
    for gene in genes:
        try:
            with open(os.path.join(out_dir, 'thoraxe_logs', gene + '.log'), encoding='utf-8',
                      errors='replace') as stream:
                lines = stream.read().splitlines()
        except OSError:
            codes['not started'] = codes.get('not started', 0) + 1
            continue
        last = lines[-1] if lines else ''
        key = last if last.startswith('exit=') else 'running or interrupted'
        codes[key] = codes.get(key, 0) + 1
        if key.startswith('exit=') and key != 'exit=0':
            failed.append(gene)
        for line in lines[-3:]:
            if line.startswith('resources='):
                wall += json.loads(line.split('=', 1)[1])['wall_s']
    print(f'{len(genes)} genes: ' + ', '.join(f'{k}: {v}' for k, v in sorted(codes.items())))
    print(f'ThorAxe time of the finished genes: {wall / 3600:.2f} CPU hours ({wall:.0f} s)')
    if failed:
        print('failed: ' + ' '.join(failed[:50]) + (' ...' if len(failed) > 50 else ''))
    with open(os.path.join(out_dir, INFO), encoding='utf-8') as stream:
        jobs = [s['job_id'] for s in json.load(stream)['submissions']]
    if jobs:
        queue = subprocess.run(['squeue', '-h', '-j', ','.join(jobs), '-o', '%T'],
                               stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                               universal_newlines=True, check=False).stdout.split()
        states = {state: queue.count(state) for state in sorted(set(queue))}
        print(f'jobs {",".join(jobs)} in queue: {states or "none"}')


def main():  # pylint: disable=too-many-statements,too-many-branches
    "Command line."
    argv = sys.argv[1:]
    thoraxe_args = []
    if '--' in argv:
        thoraxe_args = argv[argv.index('--') + 1:]
        argv = argv[:argv.index('--')]
    # "--sbatch --qos=x": argparse would take the value for an option
    for i in range(len(argv) - 1):
        if argv[i] == '--sbatch':
            argv[i:i + 2] = ['--sbatch=' + argv[i + 1], '']
    argv = [arg for arg in argv if arg != '']
    parser = argparse.ArgumentParser(
        prog='thoraxe_array.py',
        usage='%(prog)s [options] IN_DIR OUT_DIR -- [thoraxe options]\n'
              '       %(prog)s --resume OUT_DIR [--sbatch OPTIONS]\n'
              '       %(prog)s --status OUT_DIR',
        description=__doc__.split('\n\n', 2)[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dirs', nargs='+', metavar='DIR', help='IN_DIR OUT_DIR, or OUT_DIR')
    parser.add_argument('--python', default=sys.executable,
                        help='interpreter where thoraxe is installed (default: this one)')
    parser.add_argument('--seed', type=int, default=0, help='PYTHONHASHSEED (default: 0)')
    parser.add_argument('--genes', help='file with the genes to run (default: every gene '
                        'folder of IN_DIR with Ensembl/exonstable.tsv)')
    parser.add_argument('--provenance', help='JSON file describing how IN_DIR was built, '
                        'stored in thoraxe_info.json')
    parser.add_argument('--sbatch', default='', help='sbatch options, in one string')
    parser.add_argument('--dry-run', action='store_true',
                        help='check everything and print the sbatch command; write nothing')
    parser.add_argument('--resume', action='store_true', help='submit OUT_DIR again')
    parser.add_argument('--status', action='store_true', help='print the progress of OUT_DIR')
    opts = parser.parse_args(argv)

    if opts.resume or opts.status:
        if len(opts.dirs) != 1 or thoraxe_args:
            parser.error('--resume and --status take OUT_DIR only')
        out_dir = os.path.abspath(opts.dirs[0])
        if not os.path.isfile(os.path.join(out_dir, INFO)):
            raise SystemExit(f'{out_dir} was not prepared by thoraxe_array.py')
        if opts.status:
            status(out_dir)
        else:
            submit(out_dir, opts.sbatch, opts.dry_run)
        return

    if len(opts.dirs) != 2:
        parser.error('IN_DIR and OUT_DIR are needed')
    in_dir, out_dir = (os.path.abspath(d) for d in opts.dirs)
    if not os.path.isdir(in_dir):
        raise SystemExit(f'no folder {in_dir}')
    if os.path.exists(os.path.join(out_dir, INFO)):
        raise SystemExit(f'{out_dir} is already prepared: use --resume (or another OUT_DIR)')
    if os.path.isdir(out_dir) and os.listdir(out_dir):
        raise SystemExit(f'{out_dir} is not empty')

    if opts.dry_run:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = os.path.join(tmp, 'run')
            os.makedirs(run_dir)
            info = prepare(opts, in_dir, out_dir, thoraxe_args, run_dir)
            print(json.dumps({k: info[k] for k in ('input', 'thoraxe', 'species_list',
                                                   'aligner', 'python')}, indent=1))
            with open(os.path.join(run_dir, 'genes.txt'), encoding='utf-8') as stream:
                first = [next(stream).strip() for _ in range(min(5, info['input']['genes_run']))]
            print(f'first genes (largest exon tables): {" ".join(first)}')
        submit(out_dir, opts.sbatch, dry_run=True)
        return

    run_dir = os.path.join(out_dir, 'run')
    os.makedirs(run_dir)
    os.makedirs(os.path.join(out_dir, 'thoraxe_logs'))
    info = prepare(opts, in_dir, out_dir, thoraxe_args, run_dir)
    _write_info(out_dir, info)
    print(f'prepared {out_dir}: {info["input"]["genes_run"]} genes '
          f'({info["input"]["genes_without_exon_table"]} gene folders without exon table)')
    submit(out_dir, opts.sbatch)
    print(f'progress: python {os.path.abspath(__file__)} --status {out_dir}\n'
          f'resume:   python {os.path.abspath(__file__)} --resume {out_dir}')


if __name__ == '__main__':
    main()
