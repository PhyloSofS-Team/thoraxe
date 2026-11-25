"""
index: build the indexes that transcript_query --bulk-dir reads.

Each index is written to a temporary file and renamed when complete, and
existing indexes are kept: running the command again resumes an interrupted
run. Delete an index to rebuild it.
"""

import concurrent.futures
import os
import sqlite3
import time

from thoraxe.bulk_ensembl import compara, exons, sequences, tsl


def _log(message):
    print(f'[index] {message}', flush=True)


def _indexed_tsv_species(homologies_db):
    "Species whose TSV files the homology index was built from."
    con = sqlite3.connect(f'file:{homologies_db}?mode=ro', uri=True)
    try:
        row = con.execute("SELECT value FROM meta WHERE key = 'tsv_species'").fetchone()
    except sqlite3.OperationalError:
        row = None
    finally:
        con.close()
    return row[0].split(',') if row and row[0] else []


def index_compara(bulk):
    """
    Build the gene tree and homology indexes of the bulk folder.

    The homology index is rebuilt when the downloaded homology files change.
    Raise if a species that shares gene trees with the query species has no
    homology in the index (its homology file has to be downloaded).
    """
    if not os.path.isfile(bulk.trees_db()):
        start = time.time()
        compara.build_tree_index(bulk.compara_emf(), bulk.trees_db())
        _log(f'gene trees: {time.time() - start:.0f} s')

    query = bulk.query_species
    others = [species for species in bulk.compara_tsv_species() if species != query]
    db_path = bulk.homologies_db()
    if (os.path.isfile(db_path) and
            (compara.index_query_species(db_path) != query or
             _indexed_tsv_species(db_path) != others)):
        os.remove(db_path)
    if not os.path.isfile(db_path):
        start = time.time()
        compara.build_homology_index(query, bulk.compara_tsv(query),
                                     [bulk.compara_tsv(species) for species in others],
                                     db_path)
        con = sqlite3.connect(db_path)
        con.execute("INSERT INTO meta VALUES ('tsv_species', ?)", (','.join(others),))
        con.commit()
        con.close()
        _log(f'homologies of {query} (+ {len(others)} reverse files): '
             f'{time.time() - start:.0f} s')

    missing = compara.missing_homology_species(bulk.trees_db(), db_path, query)
    if missing:
        raise RuntimeError(
            f'{len(missing)} species share gene trees with {query} but have no '
            f'homology in the index: {", ".join(missing)}. Their pairs are only '
            'in their own Compara homology file: download it with "ensembl_bulk '
            f'download {bulk.root} --homology-species {" ".join(missing)} ..." '
            'and run "ensembl_bulk index" again.')


def index_species(bulk, species):
    "Build the tsl tables, exon table, genome and GTF indexes of species."
    start = time.time()
    built = []
    if bulk.missing_indexes(species):
        if not os.path.isfile(os.path.join(bulk.tsl_tables_dir(species),
                                           'transcripts.tsv')):
            tsl.build_tables(bulk.core_dir(species), bulk.tsl_tables_dir(species))
            built.append('tsl tables')
        if not (os.path.isfile(bulk.exons_db(species)) or
                os.path.isfile(bulk.exons_no_dataset(species))):
            exons.build_index(bulk.exons_raw(species), bulk.exons_db(species))
            built.append('exon table')
        fa_path = bulk.genome_fa(species)
        if not os.path.isfile(fa_path + '.fai'):
            if not os.path.isfile(fa_path):
                sequences.decompress(bulk.genome_gz(species), fa_path)
            sequences.build_fai(fa_path)
            built.append('genome')
        if not os.path.isfile(bulk.gtf_db(species)):
            sequences.build_gtf_index(bulk.gtf_gz(species), bulk.gtf_db(species))
            built.append('gtf')
    return species, built, time.time() - start


def index(bulk, species_list=None, workers=1):
    """
    Build every index of the bulk folder.

    species_list defaults to the species of bulk_info.json. workers species are
    indexed in parallel (each one needs a few GB of memory).
    """
    index_compara(bulk)
    species_list = species_list or bulk.species or []
    if bulk.query_species not in species_list:
        species_list = [bulk.query_species] + list(species_list)
    failed = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = {pool.submit(index_species, bulk, species): species
                   for species in species_list}
        for future in concurrent.futures.as_completed(futures):
            try:
                species, built, seconds = future.result()
            except Exception as err:  # pylint: disable=broad-except
                failed.append(f'{futures[future]}: {err!r}')
                continue
            _log(f'{species}: {", ".join(built) if built else "ready"} '
                 f'({seconds:.0f} s)')
    if failed:
        raise RuntimeError('indexing failed for:\n' + '\n'.join(failed))
    _log(f'{len(species_list)} species ready in {bulk.root}')
