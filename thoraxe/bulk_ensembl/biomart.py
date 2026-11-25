"""
biomart: download one exon table per species from the BioMart of a release.

transcript_query asks BioMart for the exon table of each gene. Here the same
query (same 14 attributes, same dataset) is sent once per species without the
gene filter, to the BioMart of the requested release:

* the archive site of the release (``ARCHIVE_HOSTS``) or ``--biomart-host``;
* a host is used only if its mart registry advertises "Ensembl Genes <R>", so
  that a dump never silently comes from another release;
* the answer must end with BioMart's completion stamp ``[success]``, so that a
  truncated answer is never stored.
"""

import gzip
import os
import re
import time

import requests

from thoraxe.transcript_query.transcript_query import _species2ensembldataset

BIOMART_PATH = '/biomart/martservice'
MART = 'ENSEMBL_MART_ENSEMBL'

#: Archive site of each Ensembl release (https only).
ARCHIVE_HOSTS = {
    112: 'https://may2024.archive.ensembl.org',
    113: 'https://oct2024.archive.ensembl.org',
    114: 'https://may2025.archive.ensembl.org',
    115: 'https://sep2025.archive.ensembl.org',
    116: 'https://jun2026.archive.ensembl.org',
}

#: The BioMart attributes of transcript_query, in its order.
ATTRIBUTES = [
    'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'strand',
    'ensembl_exon_id', 'exon_chrom_start', 'exon_chrom_end', 'rank',
    'cdna_coding_start', 'cdna_coding_end', 'genomic_coding_start',
    'genomic_coding_end', 'phase', 'end_phase'
]

#: Some mirrors refuse the default python-requests User-Agent.
HEADERS = {'User-Agent': 'thoraxe-ensembl-bulk'}
RETRY_WAITS = (15, 30, 60, 120, 240)
COMPLETION_STAMP = '[success]'


def xml_query(dataset, header=True):
    """
    Return the BioMart XML query of transcript_query, without gene filter.

    The completion stamp makes BioMart end a complete answer with [success].
    """
    attributes = ''.join(f'<Attribute name = "{name}" />' for name in ATTRIBUTES)
    return ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
            '<Query  virtualSchemaName = "default" formatter = "TSV" '
            f'header = "{int(header)}" uniqueRows = "0" count = "" '
            'datasetConfigVersion = "0.6" completionStamp = "1">'
            f'<Dataset name = "{dataset}" interface = "default">'
            f'{attributes}</Dataset></Query>')


def _get(url, params, timeout, stream=False, waits=RETRY_WAITS):
    "GET retried on connection errors, 429 and 5xx answers."
    error = None
    for wait in tuple(waits) + (None,):
        try:
            response = requests.get(url, params=params, headers=HEADERS,
                                    timeout=timeout, stream=stream)
            if response.status_code != 429 and response.status_code < 500:
                response.raise_for_status()
                return response
            error = f'HTTP {response.status_code}'
        except (requests.ConnectionError, requests.Timeout,
                requests.exceptions.ChunkedEncodingError) as err:
            error = type(err).__name__
        if wait is None:
            break
        print(f'[biomart] {error}; retrying in {wait} s', flush=True)
        time.sleep(wait)
    raise requests.HTTPError(f'BioMart request failed ({error}): {url}')


def served_release(host):
    "Ensembl release advertised by the mart registry of host (None if none)."
    try:
        response = _get(host + BIOMART_PATH, {'type': 'registry'}, 60,
                        waits=(15, 60))
    except requests.RequestException:
        return None
    match = re.search(r'displayName="Ensembl Genes (\d+)"', response.text)
    return int(match.group(1)) if match else None


def biomart_host(release, host=None):
    """
    Return a BioMart host serving release (host if given, else the archive).

    Raise ValueError if the host serves another release.
    """
    host = (host or ARCHIVE_HOSTS.get(int(release), '')).rstrip('/')
    if not host:
        raise ValueError(
            f'no known BioMart archive for Ensembl {release}: give the '
            'archive site with --biomart-host (see '
            'https://www.ensembl.org/info/website/archives/index.html)')
    served = served_release(host)
    if served != int(release):
        raise ValueError(
            f'{host} serves Ensembl {served}, not Ensembl {release}')
    return host


def gene_datasets(host):
    "Names of the gene datasets of the Ensembl mart of host."
    response = _get(host + BIOMART_PATH, {'type': 'datasets', 'mart': MART}, 600)
    names = set()
    for line in response.text.splitlines():
        parts = line.split('\t')
        if len(parts) > 2 and parts[0] == 'TableSet':
            names.add(parts[1])
    return names


def dataset_for_species(species):
    "BioMart gene dataset of species, as transcript_query names it."
    return _species2ensembldataset(species)[0]


def download_species(host, species, raw_path, no_dataset_path, datasets=None):
    """
    Download the exon table of species into raw_path (gzipped TSV).

    If the mart has no dataset for species, write no_dataset_path instead
    (transcript_query records such genes in errors.csv). Existing outputs are
    kept.
    """
    if os.path.isfile(raw_path) or os.path.isfile(no_dataset_path):
        return raw_path
    os.makedirs(os.path.dirname(raw_path), exist_ok=True)
    dataset = dataset_for_species(species)
    if datasets is not None and dataset not in datasets:
        with open(no_dataset_path, 'w', encoding='utf-8') as stream:
            stream.write(f'{dataset} is not a dataset of {MART} at {host}\n')
        return no_dataset_path
    start = time.time()
    part = raw_path + '.download'
    response = _get(host + BIOMART_PATH, {'query': xml_query(dataset)}, 7200,
                    stream=True)
    with open(part, 'wb') as stream:
        for chunk in response.iter_content(chunk_size=1 << 20):
            stream.write(chunk)
    _store_checked(part, raw_path, f'{species} ({dataset})')
    os.remove(part)
    print(f'[biomart] {species}: {os.path.getsize(raw_path)} bytes, '
          f'{time.time() - start:.0f} s', flush=True)
    return raw_path


def _store_checked(part, raw_path, what):
    "Check the BioMart answer in part and gzip it without the stamp."
    size = os.path.getsize(part)
    with open(part, 'rb') as stream:
        head = stream.read(200)
        stream.seek(max(0, size - 100))
        tail = stream.read()
    if head.lstrip().startswith(b'<') or head.startswith(b'Query ERROR'):
        raise RuntimeError(f'{what}: BioMart error: {head!r}')
    if not tail.rstrip(b'\n').endswith(COMPLETION_STAMP.encode()):
        raise RuntimeError(f'{what}: incomplete BioMart answer (no '
                           f'{COMPLETION_STAMP} at the end): ...{tail!r}')
    stamp_at = size - len(tail) + tail.rfind(COMPLETION_STAMP.encode())
    tmp = raw_path + '.part'
    with open(part, 'rb') as source, gzip.open(tmp, 'wb') as target:
        left = stamp_at
        while left > 0:
            chunk = source.read(min(left, 1 << 24))
            if not chunk:
                raise RuntimeError(f'{what}: short read while storing the answer')
            target.write(chunk)
            left -= len(chunk)
    os.replace(tmp, raw_path)
