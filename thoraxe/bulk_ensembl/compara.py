"""
compara: gene trees and orthologs from the Ensembl Compara bulk files.

transcript_query asks, per gene,

* ``tree.nh``: REST ``/genetree/member/id/<species>/<gene>?nh_format=full``;
* the orthologs: REST ``/homology/id/<species>/<gene>``.

Here both come from the ``protein_default`` collection of Compara:

* ``Compara.<R>.protein_default.nh.emf.gz``: every gene tree (Newick) and its
  members, indexed in ``trees.sqlite``;
* ``Compara.<R>.protein_default.homologies.tsv.gz``: the homology pairs of one
  species, indexed in ``homologies.sqlite`` for the query species.

The homology file of a species does not hold all of its pairs: in release 115,
the human file lacks every pair whose other gene belongs to
cricetulus_griseus_picr, equus_caballus, mus_musculus, ovis_aries or sus_scrofa.
Those pairs are only in the files of these species, with the human gene as
homology partner. The index therefore reads the query species file and, in
reverse, the files of the other downloaded species. `missing_homology_species`
finds the species whose file is needed for a complete index.

A tree from the EMF dump has the same leaves, labels and branch lengths as the
REST tree; the order of sibling subtrees can differ, as it does between two
REST calls.
"""

import gzip
import os
import sqlite3

from thoraxe.transcript_query.transcript_query import FORBIDDEN_SPECIES

#: Species whose homology file carries human pairs absent from the human file
#: (Ensembl release 115, see the module documentation).
DEFAULT_HOMOLOGY_SPECIES = {
    'homo_sapiens': ['cricetulus_griseus_picr', 'equus_caballus',
                     'mus_musculus', 'ovis_aries', 'sus_scrofa']
}

_BATCH = 200000


def _tmp_connection(db_path):
    "Return a connection to db_path + '.part', removing any old one."
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    tmp = db_path + '.part'
    if os.path.exists(tmp):
        os.remove(tmp)
    return sqlite3.connect(tmp), tmp


def build_tree_index(emf_path, db_path):
    """
    Index the EMF Newick dump into db_path.

    Tables: members(gene_id, species, peptide_id, gene_name, tree_id) and
    trees(tree_id, newick).
    """
    con, tmp = _tmp_connection(db_path)
    con.executescript("""
        CREATE TABLE members(gene_id TEXT, species TEXT, peptide_id TEXT,
                             gene_name TEXT, tree_id INTEGER);
        CREATE TABLE trees(tree_id INTEGER PRIMARY KEY, newick TEXT NOT NULL);
    """)
    members, trees = [], []
    tree_id = -1
    in_block = False
    state = None
    with gzip.open(emf_path, 'rt', encoding='utf-8') as stream:
        for line in stream:
            line = line.rstrip('\n')
            if line.startswith('SEQ '):
                if not in_block:
                    tree_id += 1
                    in_block = True
                fields = line.split(' ')
                # SEQ species peptide contig start end strand gene_id gene_name
                members.append((fields[7], fields[1], fields[2],
                                fields[8] if len(fields) > 8 else None, tree_id))
                state = 'seq'
            elif line == 'DATA':
                state = 'data'
            elif state == 'data':
                trees.append((tree_id, line))
                state = 'newick'
            elif line == '//':
                in_block = False
                state = None
                if len(members) > _BATCH:
                    con.executemany('INSERT INTO members VALUES (?,?,?,?,?)', members)
                    members = []
                if len(trees) > 2000:
                    con.executemany('INSERT INTO trees VALUES (?,?)', trees)
                    trees = []
    con.executemany('INSERT INTO members VALUES (?,?,?,?,?)', members)
    con.executemany('INSERT INTO trees VALUES (?,?)', trees)
    con.execute('CREATE INDEX ix_members_gene ON members(gene_id)')
    con.execute('CREATE INDEX ix_members_tree ON members(tree_id)')
    con.commit()
    con.close()
    os.replace(tmp, db_path)
    return db_path


def build_homology_index(query_species, query_tsv, other_tsvs, db_path):  # pylint: disable=too-many-locals
    """
    Index the homologies of query_species into db_path.

    Rows of query_tsv are read forward (gene_stable_id is the query gene). Rows
    of the other files are kept only when homology_species is query_species,
    and stored with the homology partner as query gene.
    Table: homology(gene_id, species, type, target_gene_id, target_species).
    """
    con, tmp = _tmp_connection(db_path)
    con.execute('CREATE TABLE homology(gene_id TEXT, species TEXT, type TEXT, '
                'target_gene_id TEXT, target_species TEXT)')
    con.execute('CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)')
    for path in [query_tsv] + list(other_tsvs):
        rows = []
        with gzip.open(path, 'rt', encoding='utf-8') as stream:
            header = stream.readline().rstrip('\n').split('\t')
            i_gene = header.index('gene_stable_id')
            i_species = header.index('species')
            i_type = header.index('homology_type')
            i_target = header.index('homology_gene_stable_id')
            i_target_species = header.index('homology_species')
            for line in stream:
                fields = line.rstrip('\n').split('\t')
                if fields[i_species] == query_species:
                    rows.append((fields[i_gene], fields[i_species], fields[i_type],
                                 fields[i_target], fields[i_target_species]))
                elif fields[i_target_species] == query_species:
                    rows.append((fields[i_target], fields[i_target_species],
                                 fields[i_type], fields[i_gene], fields[i_species]))
                if len(rows) >= _BATCH:
                    con.executemany('INSERT INTO homology VALUES (?,?,?,?,?)', rows)
                    rows = []
        con.executemany('INSERT INTO homology VALUES (?,?,?,?,?)', rows)
    # a pair present in two files is kept once, at its first position
    con.execute('DELETE FROM homology WHERE rowid NOT IN '
                '(SELECT MIN(rowid) FROM homology GROUP BY gene_id, target_gene_id)')
    con.execute('CREATE INDEX ix_hom_gene ON homology(gene_id)')
    con.execute("INSERT INTO meta VALUES ('query_species', ?)", (query_species,))
    con.commit()
    con.close()
    os.replace(tmp, db_path)
    return db_path


def ortholog_species(query_species, query_tsv, other_tsvs):
    """
    Return the species with an ortholog of a query_species gene.

    These are the species that transcript_query can query without
    --specieslist: any orthology type, at least two name parts, no hybrid
    cattle. The homology files are read as in build_homology_index.
    """
    species = set()
    for path in [query_tsv] + list(other_tsvs):
        with gzip.open(path, 'rt', encoding='utf-8') as stream:
            header = stream.readline().rstrip('\n').split('\t')
            i_species = header.index('species')
            i_type = header.index('homology_type')
            i_target_species = header.index('homology_species')
            for line in stream:
                fields = line.rstrip('\n').split('\t')
                if not fields[i_type].startswith('ortholog_'):
                    continue
                if fields[i_species] == query_species:
                    species.add(fields[i_target_species])
                elif fields[i_target_species] == query_species:
                    species.add(fields[i_species])
    return sorted(name for name in species
                  if len(name.split('_')) >= 2 and name not in FORBIDDEN_SPECIES)


def missing_homology_species(trees_db, homologies_db, query_species):
    """
    Return the species that share a gene tree with a query_species gene but
    have no homology with it in the index (their homology file is needed).
    """
    con = sqlite3.connect(f'file:{trees_db}?mode=ro', uri=True)
    try:
        in_trees = {row[0] for row in con.execute(
            'SELECT DISTINCT species FROM members WHERE tree_id IN '
            '(SELECT tree_id FROM members WHERE species = ?)', (query_species,))}
    finally:
        con.close()
    con = sqlite3.connect(f'file:{homologies_db}?mode=ro', uri=True)
    try:
        in_homologies = {row[0] for row in con.execute(
            'SELECT DISTINCT target_species FROM homology')}
    finally:
        con.close()
    return sorted(in_trees - in_homologies - {query_species})


def index_query_species(homologies_db):
    "Return the query species the homology index was built for."
    con = sqlite3.connect(f'file:{homologies_db}?mode=ro', uri=True)
    try:
        row = con.execute(
            "SELECT value FROM meta WHERE key = 'query_species'").fetchone()
    except sqlite3.OperationalError:  # index without meta table
        row = None
    finally:
        con.close()
    return row[0] if row else None


def tree_newick(trees_db, gene_id):
    "Return the Newick text of the gene tree containing gene_id (None if none)."
    con = sqlite3.connect(f'file:{trees_db}?mode=ro', uri=True)
    try:
        row = con.execute(
            'SELECT t.newick FROM members m JOIN trees t ON t.tree_id = m.tree_id '
            'WHERE m.gene_id = ? LIMIT 1', (gene_id,)).fetchone()
    finally:
        con.close()
    return row[0] if row else None


def homologies(homologies_db, gene_id):
    """
    Return the homologies of gene_id as REST /homology/id returns them.

    Each element is a dict with 'type', 'source' and 'target' ('species' and
    'id'), the fields that transcript_query uses. Index order.
    """
    con = sqlite3.connect(f'file:{homologies_db}?mode=ro', uri=True)
    try:
        rows = con.execute(
            'SELECT species, type, target_gene_id, target_species FROM homology '
            'WHERE gene_id = ? ORDER BY rowid', (gene_id,)).fetchall()
    finally:
        con.close()
    return [{'type': kind,
             'source': {'species': species, 'id': gene_id},
             'target': {'species': target_species, 'id': target}}
            for species, kind, target, target_species in rows]
