"""
exons: exonstable.tsv from the BioMart dump of each species.

The dump of a species (biomart.download_species) holds the answer that
transcript_query would get for each of its genes, one after the other. It is
indexed by gene (first column) in ``<species>.sqlite``, keeping the dump order,
so that `exon_table` returns the text transcript_query writes for a gene.
"""

import csv
import gzip
import os
import sqlite3

from thoraxe.bulk_ensembl.biomart import ATTRIBUTES

_COLUMNS = ','.join(f'c{i}' for i in range(len(ATTRIBUTES)))


def build_index(raw_path, db_path):
    """
    Index the BioMart dump raw_path into db_path.

    Tables: exons(c0..c13, the dump columns, rowid = dump order) and
    meta(key, value) with the dump header.
    """
    tmp = db_path + '.part'
    if os.path.exists(tmp):
        os.remove(tmp)
    con = sqlite3.connect(tmp)
    con.execute(f"CREATE TABLE exons ({_COLUMNS.replace(',', ' TEXT,')} TEXT)")
    con.execute('CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)')
    insert = f"INSERT INTO exons VALUES ({','.join('?' * len(ATTRIBUTES))})"
    header = None
    batch = []
    with gzip.open(raw_path, 'rt', encoding='utf-8', newline='') as stream:
        for row in csv.reader(stream, delimiter='\t', quoting=csv.QUOTE_NONE):
            if header is None:
                header = '\t'.join(row)
                continue
            if len(row) != len(ATTRIBUTES):
                raise ValueError(f'{raw_path}: row with {len(row)} fields: {row[:3]}')
            batch.append(row)
            if len(batch) >= 50000:
                con.executemany(insert, batch)
                batch = []
    con.executemany(insert, batch)
    con.execute('CREATE INDEX ix_gene ON exons (c0)')
    con.execute("INSERT INTO meta VALUES ('header', ?)", (header,))
    con.commit()
    con.close()
    os.replace(tmp, db_path)
    return db_path


def exon_table(db_path, gene_id, header=True):
    """
    Return the BioMart answer for gene_id (None if the gene has no rows).

    The header line is included if header is True. The caller applies
    transcript_query's column renaming.
    """
    con = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    try:
        rows = con.execute(f'SELECT {_COLUMNS} FROM exons WHERE c0 = ? '
                           'ORDER BY rowid', (gene_id,)).fetchall()
        header_line = con.execute(
            "SELECT value FROM meta WHERE key = 'header'").fetchone()[0]
    finally:
        con.close()
    if not rows:
        return None
    lines = [header_line] if header else []
    lines.extend('\t'.join(row) for row in rows)
    return '\n'.join(lines) + '\n'
