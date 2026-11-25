"""
Tests of transcript_query --bulk-dir on a small synthetic bulk folder.

The folder holds one gene tree with a human gene (ENSG00000000001, GENE1) and
its mouse ortholog, the Compara, core database, BioMart, GTF and genome files
of both species, and is indexed with ensembl_bulk index before the query.
"""

import gzip
import os
import sys

import pytest

from thoraxe import transcript_query
from thoraxe.bulk_ensembl import biomart, compara, download, index
from thoraxe.bulk_ensembl.layout import BulkDir

HUMAN_CHR = ('ACGTTGCA' * 15)  # 120 bases
MOUSE_CHR = ('GGATCCTA' * 10)  # 80 bases

BIOMART_HEADER = '\t'.join([
    'Gene stable ID', 'Transcript stable ID', 'Protein stable ID', 'Strand',
    'Exon stable ID', 'Exon region start (bp)', 'Exon region end (bp)',
    'Exon rank in transcript', 'cDNA coding start', 'cDNA coding end',
    'Genomic coding start', 'Genomic coding end', 'Start phase', 'End phase'
])
HUMAN_EXONS = [
    'ENSG00000000001\tENST00000000001\tENSP00000000001\t1\tENSE00000000001\t11\t30\t1\t5\t20\t15\t30\t-1\t1',
    'ENSG00000000001\tENST00000000001\tENSP00000000001\t1\tENSE00000000002\t41\t60\t2\t21\t35\t41\t55\t1\t-1',
]
MOUSE_EXONS = [
    'ENSMUSG00000000001\tENSMUST00000000001\tENSMUSP00000000001\t-1\tENSMUSE00000000001\t5\t25\t1\t1\t21\t5\t25\t0\t0',
]

HOMOLOGY_HEADER = ('gene_stable_id\tprotein_stable_id\tspecies\tidentity\t'
                   'homology_type\thomology_gene_stable_id\t'
                   'homology_protein_stable_id\thomology_species\thomology_identity')

SCHEMA = {
    'analysis': ['analysis_id', 'created', 'logic_name'],
    'attrib_type': ['attrib_type_id', 'code', 'name'],
    'transcript_attrib': ['transcript_id', 'attrib_type_id', 'value'],
    'transcript': ['transcript_id', 'gene_id', 'analysis_id', 'seq_region_id',
                   'seq_region_start', 'seq_region_end', 'seq_region_strand',
                   'display_xref_id', 'source', 'biotype', 'stable_id', 'version',
                   'is_current'],
    'xref': ['xref_id', 'external_db_id', 'dbprimary_acc', 'display_label'],
    'gene': ['gene_id', 'biotype', 'analysis_id', 'seq_region_id',
             'seq_region_start', 'seq_region_end', 'seq_region_strand',
             'stable_id', 'is_current'],
}


def _gz(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, 'wt', encoding='utf-8') as stream:
        stream.write(text)


def _fasta(name, assembly, sequence):
    lines = [f'>{name} dna:chromosome chromosome:{assembly}:{name}:1:'
             f'{len(sequence)}:1 REF']
    lines += [sequence[i:i + 60] for i in range(0, len(sequence), 60)]
    return '\n'.join(lines) + '\n'


def _gtf(rows):
    lines = ['#!genome-build test']
    for seq, feature, start, end, strand, attributes in rows:
        text = ' '.join(f'{key} "{value}";' for key, value in attributes.items())
        lines.append('\t'.join([seq, 'ensembl', feature, str(start), str(end), '.',
                                strand, '.', text]))
    return '\n'.join(lines) + '\n'


def _core(root, species, tables, name=None):
    name = name or f'{species}_core_115_1'
    core = os.path.join(root, 'mysql', name)
    schema = []
    for table, columns in SCHEMA.items():
        schema.append(f'CREATE TABLE `{table}` (')
        schema += [f'  `{column}` varchar(40) DEFAULT NULL,' for column in columns]
        schema += [f'  PRIMARY KEY (`{columns[0]}`)', ') ENGINE=MyISAM;']
        rows = tables.get(table, [])
        _gz(os.path.join(core, f'{table}.txt.gz'),
            ''.join('\t'.join(row) + '\n' for row in rows))
    _gz(os.path.join(core, f'{name}.sql.gz'), '\n'.join(schema) + '\n')


def make_bulk(root, with_mouse_homologies=True):
    "Write the synthetic bulk folder into root."
    _gz(os.path.join(root, 'compara', 'emf', 'Compara.115.protein_default.nh.emf.gz'),
        'SEQ homo_sapiens ENSP00000000001 1 11 60 1 ENSG00000000001 GENE1\n'
        'SEQ mus_musculus ENSMUSP00000000001 2 5 25 -1 ENSMUSG00000000001 Gene1\n'
        'DATA\n(ENSP00000000001:0.1,ENSMUSP00000000001:0.2);\n//\n')
    _gz(os.path.join(root, 'compara', 'tsv', 'homo_sapiens',
                     'Compara.115.protein_default.homologies.tsv.gz'),
        HOMOLOGY_HEADER + '\n' +
        'ENSG00000000001\tENSP00000000001\thomo_sapiens\t50\twithin_species_paralog'
        '\tENSG00000000002\tENSP00000000002\thomo_sapiens\t50\n')
    if with_mouse_homologies:  # the human-mouse pair is only in the mouse file
        _gz(os.path.join(root, 'compara', 'tsv', 'mus_musculus',
                         'Compara.115.protein_default.homologies.tsv.gz'),
            HOMOLOGY_HEADER + '\n' +
            'ENSMUSG00000000001\tENSMUSP00000000001\tmus_musculus\t90\t'
            'ortholog_one2one\tENSG00000000001\tENSP00000000001\thomo_sapiens\t90\n')

    _core(root, 'homo_sapiens', {
        'analysis': [['7', '\\N', 'ensembl'], ['8', '\\N', 'havana_lnc']],
        'attrib_type': [['1', 'TSL', 'tsl'], ['2', 'other', 'other']],
        'transcript_attrib': [['101', '1', 'tsl1 (assigned to previous version 1)'],
                              ['102', '2', 'x']],
        'transcript': [
            ['102', '2', '8', '10', '50', '90', '-1', '\\N', 'havana', 'lncRNA',
             'ENST00000000002', '1', '1'],
            ['101', '1', '7', '10', '11', '60', '1', '501', 'ensembl',
             'protein_coding', 'ENST00000000001', '2', '1'],
            ['103', '1', '7', '10', '11', '60', '1', '\\N', 'ensembl',
             'protein_coding', 'ENST00000000009', '1', '0'],
        ],
        'xref': [['501', '1', 'X', 'GENE1-201']],
        'gene': [['1', 'protein_coding', '7', '10', '11', '60', '1', 'ENSG00000000001', '1'],
                 ['2', 'lncRNA', '8', '10', '50', '90', '-1', 'ENSG00000000002', '1']],
    })
    _core(root, 'mus_musculus', {
        'analysis': [['7', '\\N', 'ensembl']],
        'attrib_type': [['1', 'TSL', 'tsl']],
        'transcript': [['201', '11', '7', '20', '5', '25', '-1', '\\N', 'ensembl',
                        'protein_coding', 'ENSMUST00000000001', '3', '1']],
        'gene': [['11', 'protein_coding', '7', '20', '5', '25', '-1',
                  'ENSMUSG00000000001', '1']],
    }, name='mus_musculus_core_62_115_1')  # Ensembl Genomes naming (EG 62 = 115)

    _gz(os.path.join(root, 'exons', 'raw', 'homo_sapiens.tsv.gz'),
        '\n'.join([BIOMART_HEADER] + HUMAN_EXONS) + '\n')
    _gz(os.path.join(root, 'exons', 'raw', 'mus_musculus.tsv.gz'),
        '\n'.join([BIOMART_HEADER] + MOUSE_EXONS) + '\n')

    _gz(os.path.join(root, 'gtf', 'homo_sapiens', 'Homo_sapiens.GRCh38.115.gtf.gz'),
        _gtf([('1', 'gene', 11, 60, '+', {'gene_id': 'ENSG00000000001',
                                         'gene_name': 'GENE1',
                                         'gene_biotype': 'protein_coding'}),
              ('1', 'exon', 11, 30, '+', {'gene_id': 'ENSG00000000001',
                                         'exon_id': 'ENSE00000000001'}),
              ('1', 'exon', 41, 60, '+', {'gene_id': 'ENSG00000000001',
                                         'exon_id': 'ENSE00000000002'}),
              ('1', 'gene', 50, 90, '-', {'gene_id': 'ENSG00000000002',
                                         'gene_biotype': 'lncRNA'}),
              ('1', 'exon', 50, 70, '-', {'gene_id': 'ENSG00000000002',
                                         'exon_id': 'ENSE00000000003'})]))
    _gz(os.path.join(root, 'gtf', 'mus_musculus', 'Mus_musculus.GRCm39.62.gtf.gz'),
        _gtf([('2', 'gene', 5, 25, '-', {'gene_id': 'ENSMUSG00000000001',
                                        'gene_name': 'Gene1'}),
              ('2', 'exon', 5, 25, '-', {'gene_id': 'ENSMUSG00000000001',
                                        'exon_id': 'ENSMUSE00000000001'})]))
    _gz(os.path.join(root, 'genome', 'homo_sapiens',
                     'Homo_sapiens.GRCh38.dna.toplevel.fa.gz'),
        _fasta('1', 'GRCh38', HUMAN_CHR))
    _gz(os.path.join(root, 'genome', 'mus_musculus',
                     'Mus_musculus.GRCm39.dna.toplevel.fa.gz'),
        _fasta('2', 'GRCm39', MOUSE_CHR))

    bulk = BulkDir(root, release=115)
    bulk.write_info(release=115, query_species='homo_sapiens',
                    species=['homo_sapiens', 'mus_musculus'])
    return bulk


def _revcomp(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(sequence))


@pytest.fixture(scope='module')
def bulk(tmp_path_factory):
    "Indexed synthetic bulk folder."
    folder = make_bulk(str(tmp_path_factory.mktemp('bulk')))
    index.index(folder)
    return folder


def _run_query(monkeypatch, cwd, *args):
    monkeypatch.chdir(cwd)
    monkeypatch.setattr(sys, 'argv', ['transcript_query'] + list(args))
    assert transcript_query.main() is None


def _read(*path):
    with open(os.path.join(*path), encoding='utf-8', newline='') as stream:
        return stream.read()


def test_bulk_folder(bulk, tmp_path, monkeypatch):
    "transcript_query --bulk-dir writes the five files transcript_query writes."
    _run_query(monkeypatch, str(tmp_path), 'ENSG00000000001', '--bulk-dir',
               bulk.root, '-o', '1:n')
    folder = os.path.join(str(tmp_path), 'ENSG00000000001', 'Ensembl')
    assert sorted(os.listdir(folder)) == [
        'ensembl_version.csv', 'exonstable.tsv', 'sequences.fasta', 'tree.nh',
        'tsl.csv']

    assert _read(folder, 'ensembl_version.csv').splitlines()[1].endswith(',62')
    assert _read(folder, 'tree.nh') == '(ENSP00000000001:0.1,ENSMUSP00000000001:0.2);'

    assert _read(folder, 'tsl.csv') == (
        'Species,Name,TranscriptID,Source,ExperimentSource,Biotype,Flags,Version\r\n'
        'homo_sapiens,GENE1-201,ENST00000000001,ensembl,ensembl,protein_coding,'
        '1 (assigned to previous version 1),2\r\n'
        'homo_sapiens,NA,ENST00000000002,havana,havana_lnc,lncRNA,NA,1\r\n'
        'mus_musculus,NA,ENSMUST00000000001,ensembl,ensembl,protein_coding,NA,3\r\n')

    table = _read(folder, 'exonstable.tsv').splitlines()
    assert table[0].split('\t')[:5] == ['GeneID', 'TranscriptID', 'ProteinID',
                                        'Strand', 'ExonID']
    assert table[1:] == HUMAN_EXONS + MOUSE_EXONS

    label = '{args.species}:{args.genename}'
    assert _read(folder, 'sequences.fasta') == ''.join([
        f'>{label} ENSE00000000001 chromosome:GRCh38:1:11:30:1\n{HUMAN_CHR[10:30]}\n',
        f'>{label} ENSE00000000002 chromosome:GRCh38:1:41:60:1\n{HUMAN_CHR[40:60]}\n',
        f'>{label} ENSE00000000003 chromosome:GRCh38:1:50:70:-1\n'
        f'{_revcomp(HUMAN_CHR[49:70])}\n',
        '>mus_musculus:ENSMUSG00000000001 ENSMUSE00000000001 '
        f'chromosome:GRCm39:2:5:25:-1\n{_revcomp(MOUSE_CHR[4:25])}\n'])


def test_gene_name_and_species_list(bulk, tmp_path, monkeypatch):
    "A gene name is resolved from the GTF; -l restricts the orthologs."
    _run_query(monkeypatch, str(tmp_path), 'GENE1', '--bulk-dir', bulk.root,
               '-o', '1:1', '-l', 'homo_sapiens,mus_musculus')
    folder = os.path.join(str(tmp_path), 'GENE1', 'Ensembl')
    assert _read(folder, 'exonstable.tsv').count('ENSMUSG00000000001') == 1

    monkeypatch.setattr(sys, 'argv', ['transcript_query', 'ENSG00000000001',
                                      '--bulk-dir', bulk.root, '-l',
                                      'homo_sapiens,rattus_norvegicus'])
    with pytest.raises(Exception, match='There are not'):
        transcript_query.main()


def test_missing_homology_species(tmp_path):
    "The index refuses a homology index that misses species of the trees."
    folder = make_bulk(str(tmp_path), with_mouse_homologies=False)
    with pytest.raises(RuntimeError, match='mus_musculus'):
        index.index_compara(folder)
    assert compara.missing_homology_species(
        folder.trees_db(), folder.homologies_db(), 'homo_sapiens') == ['mus_musculus']


def test_ortholog_species(tmp_path):
    "The default download covers the species with orthologs of the query species."
    human = str(tmp_path / 'human.tsv.gz')
    _gz(human, HOMOLOGY_HEADER + '\n' +
        'ENSG1\tP1\thomo_sapiens\t50\twithin_species_paralog\tENSG2\tP2\thomo_sapiens\t50\n'
        'ENSG1\tP1\thomo_sapiens\t50\tother_paralog\tENSRNOG1\tP3\trattus_norvegicus\t50\n'
        'ENSG1\tP1\thomo_sapiens\t90\tortholog_one2many\tENSPTRG1\tP4\tpan_troglodytes\t90\n'
        'ENSG1\tP1\thomo_sapiens\t90\tortholog_one2one\tX1\tP5\tbos_indicus_hybrid\t90\n')
    mouse = str(tmp_path / 'mouse.tsv.gz')
    _gz(mouse, HOMOLOGY_HEADER + '\n' +
        'ENSMUSG1\tP6\tmus_musculus\t90\tortholog_one2one\tENSG1\tP1\thomo_sapiens\t90\n'
        'ENSMUSG1\tP6\tmus_musculus\t90\tortholog_one2one\tENSRNOG1\tP3\trattus_norvegicus\t90\n')
    assert compara.ortholog_species('homo_sapiens', human, [mouse]) == [
        'mus_musculus', 'pan_troglodytes']


def test_homology_pair_in_two_files(tmp_path):
    "A pair found in the query species file and in a reverse file is kept once."
    folder = make_bulk(str(tmp_path))
    _gz(folder.compara_tsv('homo_sapiens'),
        HOMOLOGY_HEADER + '\n' +
        'ENSG00000000001\tENSP00000000001\thomo_sapiens\t90\tortholog_one2one\t'
        'ENSMUSG00000000001\tENSMUSP00000000001\tmus_musculus\t90\n')
    index.index_compara(folder)
    rows = compara.homologies(folder.homologies_db(), 'ENSG00000000001')
    assert [row['target']['id'] for row in rows] == ['ENSMUSG00000000001']


def test_bsd_sum(tmp_path, monkeypatch):
    "The Python BSD checksum matches sum(1)."
    path = str(tmp_path / 'data.bin')
    with open(path, 'wb') as stream:
        stream.write(bytes(range(256)) * 9)
    expected = download.bsd_sum(path)
    monkeypatch.setattr(download.shutil, 'which', lambda name: None)
    assert download.bsd_sum(path) == expected
    assert expected[1] == 3  # 2304 bytes in 1 KiB blocks


def test_biomart_completion_stamp(tmp_path):
    "BioMart answers are stored without the stamp, and refused without it."
    part = str(tmp_path / 'answer')
    raw = str(tmp_path / 'raw.tsv.gz')
    with open(part, 'w', encoding='utf-8') as stream:
        stream.write(BIOMART_HEADER + '\n' + HUMAN_EXONS[0] + '\n[success]\n')
    biomart._store_checked(part, raw, 'test')  # pylint: disable=protected-access
    with gzip.open(raw, 'rt', encoding='utf-8') as stream:
        assert stream.read() == BIOMART_HEADER + '\n' + HUMAN_EXONS[0] + '\n'

    with open(part, 'w', encoding='utf-8') as stream:
        stream.write(BIOMART_HEADER + '\n' + HUMAN_EXONS[0][:20])
    with pytest.raises(RuntimeError, match='incomplete'):
        biomart._store_checked(part, str(tmp_path / 'other.tsv.gz'), 'test')  # pylint: disable=protected-access


def test_biomart_query():
    "The dump query is transcript_query's query without the gene filter."
    query = biomart.xml_query('hsapiens_gene_ensembl')
    assert 'Filter' not in query
    assert query.count('<Attribute') == 14
    assert 'completionStamp = "1"' in query
