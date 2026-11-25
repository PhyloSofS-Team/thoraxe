"""
download: fetch the bulk files of one Ensembl release.

From the FTP site (``https://ftp.ensembl.org/pub/release-<R>/``), each file is
checked against the CHECKSUMS file of its directory (BSD ``sum``, as Ensembl
publishes it):

* Compara: the protein_default gene trees (EMF) and the homology TSV of the
  query species and of the --homology-species;
* for each species: the complete GTF, the unmasked toplevel genome and the core
  database tables that tsl.csv needs. By default, the species are the query
  species and every species with an ortholog of it in the Compara homology
  files (199 species for homo_sapiens in Ensembl 115, about 120 GB).

From BioMart: one exon table per species (see biomart).

Existing files that match their checksum are kept, so an interrupted download
can be resumed by running the same command again.
"""

import concurrent.futures
import math
import os
import re
import shutil
import subprocess
import time
import urllib.error
import urllib.request

from thoraxe.bulk_ensembl import biomart, compara
from thoraxe.bulk_ensembl.layout import eg_version
from thoraxe.bulk_ensembl.tsl import TABLES

FTP_ROOT = 'https://ftp.ensembl.org/pub'
MAX_CONNECTIONS = 4


def release_url(release):
    "Root of the FTP folder of an Ensembl release."
    return f'{FTP_ROOT}/release-{release}'


def _urlopen(url, attempts=5, timeout=120):
    "urlopen retried on network errors and 5xx answers."
    for attempt in range(attempts):
        try:
            return urllib.request.urlopen(url, timeout=timeout)  # pylint: disable=consider-using-with
        except urllib.error.HTTPError as err:
            if err.code < 500 or attempt == attempts - 1:
                raise
        except (urllib.error.URLError, OSError):
            if attempt == attempts - 1:
                raise
        time.sleep(10 * (attempt + 1))
    raise RuntimeError(f'cannot open {url}')


def read_text(url):
    "Return the text at url."
    with _urlopen(url) as response:
        return response.read().decode('utf-8')


def list_dir(url):
    "Return the names listed in an FTP folder (HTML index) at url."
    return re.findall(r'href="([^"?/][^"]*?)/?"', read_text(url.rstrip('/') + '/'))


def parse_checksums(text):
    """
    Return {file name: (checksum, blocks)} from a CHECKSUMS file.

    >>> parse_checksums('45573 45635 a.gz')
    {'a.gz': (45573, 45635)}
    """
    checksums = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 3:
            checksums[parts[2]] = (int(parts[0]), int(parts[1]))
    return checksums


def bsd_sum(path):
    "BSD sum(1) of path: 16-bit rotating checksum and size in 1 KiB blocks."
    if shutil.which('sum'):
        result = subprocess.run(['sum', path], stdout=subprocess.PIPE,
                                universal_newlines=True, check=True)
        checksum, blocks = result.stdout.split()[:2]
        return int(checksum), int(blocks)
    checksum, size = 0, 0
    with open(path, 'rb') as stream:
        while True:
            chunk = stream.read(1 << 20)
            if not chunk:
                break
            size += len(chunk)
            for byte in chunk:
                checksum = ((checksum >> 1) + ((checksum & 1) << 15) + byte) & 0xFFFF
    return checksum, math.ceil(size / 1024)


def fetch(url, dest, expected):
    """
    Download url into dest unless dest already has the expected checksum.

    Return 'kept' or 'downloaded'. Raise if the download does not match.
    """
    if os.path.isfile(dest) and bsd_sum(dest) == expected:
        return 'kept'
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    part = dest + '.part'
    for attempt in range(3):
        with _urlopen(url, timeout=600) as response, open(part, 'wb') as stream:
            shutil.copyfileobj(response, stream, 1 << 20)
        if bsd_sum(part) == expected:
            os.replace(part, dest)
            return 'downloaded'
        print(f'[ftp] checksum mismatch for {url} (attempt {attempt + 1})', flush=True)
    raise RuntimeError(f'checksum mismatch for {url}')


class FtpPlan:
    "The FTP files of a bulk folder: (url, destination, expected checksum)."

    def __init__(self, bulk):
        self.bulk = bulk
        self.base = release_url(bulk.release)
        self.jobs = []
        self._checksums = {}

    def checksums(self, folder):
        "CHECKSUMS of an FTP folder (relative to the release), kept locally too."
        if folder not in self._checksums:
            self._checksums[folder] = read_text(f'{self.base}/{folder}/CHECKSUMS')
        return parse_checksums(self._checksums[folder])

    def add(self, folder, name, dest):
        "Add folder/name, saved as dest, with its CHECKSUMS next to it."
        expected = self.checksums(folder).get(name)
        if expected is None:
            raise ValueError(f'{name} is not in {self.base}/{folder}/CHECKSUMS')
        self.jobs.append((f'{self.base}/{folder}/{name}', dest, expected))
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        with open(os.path.join(os.path.dirname(dest), 'CHECKSUMS'), 'w',
                  encoding='utf-8') as stream:
            stream.write(self._checksums[folder])

    def add_compara(self, species_list):
        "Compara gene trees and the homology TSV of each species."
        release = self.bulk.release
        self.add('emf/ensembl-compara/homologies',
                 f'Compara.{release}.protein_default.nh.emf.gz',
                 self.bulk.compara_emf())
        for species in species_list:
            self.add(f'tsv/ensembl-compara/homologies/{species}',
                     f'Compara.{release}.protein_default.homologies.tsv.gz',
                     self.bulk.compara_tsv(species))

    def add_species(self, species, core_dirs):
        "GTF, genome and core database tables of species."
        release = self.bulk.release
        folder = f'gtf/{species}'
        names = [name for name in list_dir(f'{self.base}/{folder}')
                 if re.fullmatch(rf'[^/]+\.({release}|{eg_version(release)})\.gtf\.gz',
                                 name)]
        if len(names) != 1:
            raise ValueError(f'expected one complete GTF in {folder}, found {names}')
        self.add(folder, names[0], os.path.join(self.bulk.root, 'gtf', species, names[0]))

        folder = f'fasta/{species}/dna'
        names = [name for name in list_dir(f'{self.base}/{folder}')
                 if name.endswith('.dna.toplevel.fa.gz')]
        if len(names) != 1:
            raise ValueError(f'expected one toplevel genome in {folder}, found {names}')
        self.add(folder, names[0],
                 os.path.join(self.bulk.root, 'genome', species, names[0]))

        core = core_dirs.get(species)
        if core is None:
            raise ValueError(f'no core database dump of {species} in {self.base}/mysql')
        for name in [f'{table}.txt.gz' for table in TABLES] + [f'{core}.sql.gz']:
            self.add(f'mysql/{core}', name,
                     os.path.join(self.bulk.mysql_dir(), core, name))

    def core_dirs(self):
        "Core database dump folder of each species of the release."
        release = self.bulk.release
        cores = {}
        for name in list_dir(f'{self.base}/mysql'):
            match = re.fullmatch(
                rf'([a-z0-9_]+)_core_(?:{eg_version(release)}_)?{release}_\d+', name)
            if match:
                cores[match.group(1)] = name
        return cores

    def run(self, workers=MAX_CONNECTIONS):
        "Download every planned file; return {status: count}."
        counts = {}
        failed = []
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=max(1, min(workers, MAX_CONNECTIONS))) as pool:
            futures = {pool.submit(fetch, url, dest, expected): url
                       for url, dest, expected in self.jobs}
            for future in concurrent.futures.as_completed(futures):
                try:
                    status = future.result()
                except Exception as err:  # pylint: disable=broad-except
                    failed.append(f'{futures[future]}: {err}')
                    continue
                counts[status] = counts.get(status, 0) + 1
                print(f'[ftp] {status} {futures[future]}', flush=True)
        if failed:
            raise RuntimeError('failed downloads:\n' + '\n'.join(failed))
        return counts


def download(bulk, species_list, homology_species, biomart_host=None,  # pylint: disable=too-many-locals
             workers=MAX_CONNECTIONS):
    """
    Download the bulk files of species_list into the bulk folder.

    species_list None means every species with an ortholog of the query
    species in Compara, the species transcript_query uses without
    --specieslist; the query species is always added. homology_species are the
    species whose Compara homology TSV is read in reverse (see compara).
    """
    start = time.time()
    query = bulk.query_species
    reverse = [species for species in homology_species if species != query]
    plan = FtpPlan(bulk)
    plan.add_compara([query] + reverse)
    if species_list is None:  # the Compara files tell which species to download
        plan.run(workers)
        species_list = compara.ortholog_species(
            query, bulk.compara_tsv(query),
            [bulk.compara_tsv(species) for species in reverse])
        print(f'[bulk] {len(species_list)} species have orthologs of {query} genes',
              flush=True)
        plan = FtpPlan(bulk)
    species_list = sorted(set(species_list) | {query})
    core_dirs = plan.core_dirs()
    for species in species_list:
        plan.add_species(species, core_dirs)
    counts = plan.run(workers)
    print(f'[ftp] {counts} in {time.time() - start:.0f} s', flush=True)

    todo = [species for species in species_list
            if not (os.path.isfile(bulk.exons_raw(species)) or
                    os.path.isfile(bulk.exons_no_dataset(species)))]
    host = None
    if todo:
        host = biomart.biomart_host(bulk.release, biomart_host)
        datasets = biomart.gene_datasets(host)
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=max(1, min(workers, MAX_CONNECTIONS))) as pool:
            futures = [pool.submit(biomart.download_species, host, species,
                                   bulk.exons_raw(species),
                                   bulk.exons_no_dataset(species), datasets)
                       for species in todo]
            for future in futures:
                future.result()
    no_dataset = [species for species in species_list
                  if os.path.isfile(bulk.exons_no_dataset(species))]
    if no_dataset:
        print(f'[biomart] no gene dataset for: {", ".join(no_dataset)} (their '
              'genes will be recorded in errors.csv)', flush=True)
    info = bulk.read_info()
    bulk.write_info(release=bulk.release, query_species=bulk.query_species,
                    species=sorted(set(species_list) | set(info.get('species', []))),
                    homology_species=sorted(set(homology_species) |
                                            set(info.get('homology_species', []))),
                    ftp=release_url(bulk.release),
                    biomart_host=host or info.get('biomart_host'))
    print(f'[bulk] download done in {time.time() - start:.0f} s', flush=True)
