import collections
import os
import shutil
import sys
import subprocess
import warnings
import pytest

import pandas as pd
import thoraxe
from thoraxe import utils
from thoraxe import add_transcripts
from thoraxe.transcript_info import read_exon_file

Args = collections.namedtuple('Args', ['input', 'ensembl'])


# https://stackoverflow.com/a/12514470/1846377
def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


@pytest.fixture(scope='module')
def copy_mapk8_dir():
    copied_dir = 'MAPK8_copy'
    yield copied_dir
    os.path.isdir(copied_dir) and shutil.rmtree(copied_dir)


@pytest.fixture(scope='module')
def set_out_dir():
    out_dir = 'tmp'
    yield out_dir
    os.path.isdir(out_dir) and shutil.rmtree(out_dir)


def test_add_transcripts(monkeypatch, request, set_out_dir, copy_mapk8_dir):
    # Paths
    in_dir = 'MAPK8_copy'
    out_dir = 'tmp'
    filename = request.module.__file__
    data_path = os.path.join(os.path.dirname(filename), 'data')
    original_dir = os.path.join(data_path, 'MAPK8')
    copytree(original_dir, in_dir)
    args = Args(os.path.join(data_path, 'user_transcript.csv'),
                os.path.join(in_dir, "Ensembl"))
    new_paths = add_transcripts.get_output_paths(args)
    original_paths = add_transcripts.get_output_paths(
        Args('', os.path.join(original_dir, "Ensembl")))

    # Add the user transcripts
    monkeypatch.setattr(sys, 'argv',
                        ['add_transcripts', args.input, args.ensembl])
    assert add_transcripts.main() is None  # Run it

    # Test that data has been add
    assert (len(add_transcripts.read_fasta(new_paths.sequences)) == len(
        add_transcripts.read_fasta(original_paths.sequences)) + 10)
    assert (len(read_exon_file(new_paths.exontable).index) == len(
        read_exon_file(original_paths.exontable).index) + 10)
    assert (len(pd.read_csv(
        new_paths.tsl).index) == len(pd.read_csv(original_paths.tsl).index) +
            1)

    # Set aligner to run ThorAxe
    aligner = 'ProGraphMSA'
    if utils.windows.is_windows():
        print(f"WSL: {shutil.which('wsl')}")
        if shutil.which('wsl') is not None:  # ProGraphSA in the WSL
            status, _ = subprocess.getstatusoutput(
                ['wsl', 'ProGraphMSA', '-h'])
            if status == 0:
                aligner = 'wsl ProGraphMSA'
            else:
                warnings.warn(
                    'You should install ProGraphMSA in the Windows '
                    'Subsystem for Linux to run ThorAxe on the tests.')
        else:  # AppVeyor
            aligner = ""

    # Run ThorAxe with user-defined transcripts
    if aligner:
        monkeypatch.setattr(sys, 'argv', [
            'thoraxe', '-i', in_dir, '-o', out_dir, '-y', '--plot_chimerics',
            '-a', aligner
        ])
        assert not os.path.isdir('tmp')
        assert thoraxe.main() is None  # Run it
        assert os.path.isdir('tmp')
        df = pd.read_csv(os.path.join('tmp', 'thoraxe', 's_exon_table.csv'))
        assert any(df.TranscriptIDCluster == 'ENST00000360332/new_transcript')