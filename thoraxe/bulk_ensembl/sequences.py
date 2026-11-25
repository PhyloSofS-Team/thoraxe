"""
sequences: sequences.fasta from the genome FASTA and the GTF of each species.

transcript_query writes one record per exon returned by REST
``/overlap/id/<gene>?feature=exon`` (every exon, of any gene and biotype, that
overlaps the gene), with the sequence from REST ``/sequence/id`` and the
description ``<coord_system>:<assembly>:<seq_region>:<start>:<end>:<strand>``.

Here the exons and the gene spans come from the GTF (indexed in SQLite) and the
sequences from the unmasked toplevel genome, decompressed once and indexed with
a samtools-compatible ``.fai`` file. The sequence of an exon is the genomic
slice of its span, reverse complemented on the minus strand (whole exon,
UTRs included). The coordinate system and assembly of each sequence region come
from the genome FASTA headers.

Bulk files hold every exon, whereas REST ``/sequence/id`` silently skips some
exon ids (e.g. of some lncRNA genes overlapping the queried gene). The extra
records do not change thoraxe's results: thoraxe only reads the exons listed in
exonstable.tsv.
"""

import gzip
import os
import re
import sqlite3

_COMPLEMENT = str.maketrans('ACGTNacgtnRYKMSWBDHVrykmswbdhv',
                            'TGCANtgcanYRMKSWVHDBYRMKSWVHDB')
_GTF_ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')


def reverse_complement(sequence):
    """
    Return the reverse complement of a DNA sequence.

    >>> reverse_complement('AACGTN')
    'NACGTT'
    """
    return sequence.translate(_COMPLEMENT)[::-1]


# ----------------------------------------------------------------- genome


def decompress(gz_path, fa_path):
    "Decompress gz_path into fa_path."
    os.makedirs(os.path.dirname(fa_path), exist_ok=True)
    tmp = fa_path + '.part'
    with gzip.open(gz_path, 'rb') as source, open(tmp, 'wb') as target:
        while True:
            chunk = source.read(1 << 24)
            if not chunk:
                break
            target.write(chunk)
    os.replace(tmp, fa_path)
    return fa_path


def build_fai(fa_path):  # pylint: disable=too-many-locals
    """
    Index fa_path: write <fa>.regions and <fa>.fai.

    ``.fai`` (name, length, offset, line bases, line width) is the samtools
    index. ``.regions`` keeps ``coord_system:assembly`` of each sequence, from
    the third word of its FASTA header. ``.fai`` is written last: it marks the
    genome as indexed.
    """
    rows, regions = [], []
    with open(fa_path, 'rb') as stream:
        name = None
        offset = line_bases = line_width = length = 0
        cs_assembly = ''
        while True:
            line = stream.readline()
            if not line:
                break
            if line.startswith(b'>'):
                if name is not None:
                    rows.append((name, length, offset, line_bases, line_width))
                    regions.append((name, cs_assembly))
                words = line[1:].split()
                name = words[0].decode()
                fields = words[2].decode().split(':') if len(words) > 2 else []
                cs_assembly = ':'.join(fields[:2]) if len(fields) >= 2 else ''
                offset = stream.tell()
                length = line_bases = line_width = 0
                continue
            stripped = line.rstrip(b'\r\n')
            if line_bases == 0 and stripped:
                line_bases = len(stripped)
                line_width = len(line)
            length += len(stripped)
        if name is not None:
            rows.append((name, length, offset, line_bases, line_width))
            regions.append((name, cs_assembly))
    with open(fa_path + '.regions.part', 'w', encoding='utf-8') as out:
        for name, cs_assembly in regions:
            out.write(f'{name}\t{cs_assembly}\n')
    os.replace(fa_path + '.regions.part', fa_path + '.regions')
    with open(fa_path + '.fai.part', 'w', encoding='utf-8') as out:
        for row in rows:
            out.write('\t'.join(str(value) for value in row) + '\n')
    os.replace(fa_path + '.fai.part', fa_path + '.fai')
    return fa_path + '.fai'


class Genome:
    "Random access to an indexed genome FASTA."

    def __init__(self, fa_path):
        self.index = {}
        with open(fa_path + '.fai', encoding='utf-8') as stream:
            for line in stream:
                name, length, offset, line_bases, line_width = line.split('\t')
                self.index[name] = (int(length), int(offset), int(line_bases),
                                    int(line_width))
        self.cs_assembly = {}
        with open(fa_path + '.regions', encoding='utf-8') as stream:
            for line in stream:
                name, cs_assembly = line.rstrip('\n').split('\t')
                self.cs_assembly[name] = cs_assembly
        self._stream = open(fa_path, 'rb')  # pylint: disable=consider-using-with

    def close(self):
        "Close the FASTA file."
        self._stream.close()

    def fetch(self, name, start, end):
        "Return the forward strand sequence of name:start-end (1-based, inclusive)."
        length, offset, line_bases, line_width = self.index[name]
        if not 1 <= start <= end <= length:
            raise ValueError(f'{name}:{start}-{end} is outside 1..{length}')
        first = offset + ((start - 1) // line_bases) * line_width + (start - 1) % line_bases
        last = offset + ((end - 1) // line_bases) * line_width + (end - 1) % line_bases
        self._stream.seek(first)
        raw = self._stream.read(last - first + 1)
        return raw.replace(b'\n', b'').replace(b'\r', b'').decode('ascii')


# ----------------------------------------------------------------- GTF


def build_gtf_index(gtf_path, db_path):
    """
    Index the genes and exons of a GTF file into db_path.

    Tables: gene(gene_id, seq, start, end, strand, biotype, name) and
    exon(exon_id, seq, start, end, strand), one row per exon id.
    """
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    tmp = db_path + '.part'
    if os.path.exists(tmp):
        os.remove(tmp)
    con = sqlite3.connect(tmp)
    con.executescript("""
        CREATE TABLE gene(gene_id TEXT PRIMARY KEY, seq TEXT, start INT, end INT,
                          strand INT, biotype TEXT, name TEXT);
        CREATE TABLE exon(exon_id TEXT PRIMARY KEY, seq TEXT, start INT, end INT,
                          strand INT);
    """)
    genes, exons = [], {}
    with gzip.open(gtf_path, 'rt', encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] not in ('gene', 'exon'):
                continue
            attributes = dict(_GTF_ATTRIBUTE.findall(fields[8]))
            strand = 1 if fields[6] == '+' else -1
            if fields[2] == 'gene':
                genes.append((attributes['gene_id'], fields[0], int(fields[3]),
                               int(fields[4]), strand,
                               attributes.get('gene_biotype'),
                               attributes.get('gene_name')))
            elif attributes['exon_id'] not in exons:
                exons[attributes['exon_id']] = (attributes['exon_id'], fields[0],
                                                int(fields[3]), int(fields[4]),
                                                strand)
    con.executemany('INSERT INTO gene VALUES (?,?,?,?,?,?,?)', genes)
    con.executemany('INSERT INTO exon VALUES (?,?,?,?,?)', exons.values())
    con.execute('CREATE INDEX exon_span ON exon(seq, start, end)')
    con.execute('CREATE INDEX gene_span ON gene(seq, start, end)')
    con.execute('CREATE INDEX gene_name ON gene(name)')
    con.commit()
    con.close()
    os.replace(tmp, db_path)
    return db_path


def genes_named(db_path, name):
    "Return the ids of the genes called name in the GTF (case-insensitive)."
    con = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    try:
        exact = [row[0] for row in con.execute(
            'SELECT gene_id FROM gene WHERE name = ? ORDER BY gene_id', (name,))]
        if exact:
            return exact
        return [row[0] for row in con.execute(
            'SELECT gene_id FROM gene WHERE upper(name) = upper(?) '
            'ORDER BY gene_id', (name,))]
    finally:
        con.close()


class SpeciesSequences:
    "Genome and GTF index of one species, opened once for many genes."

    def __init__(self, fa_path, gtf_db):
        self.genome = Genome(fa_path)
        self._con = sqlite3.connect(f'file:{gtf_db}?mode=ro', uri=True)

    def close(self):
        "Close the genome and the GTF index."
        self.genome.close()
        self._con.close()

    def exon_records(self, gene_id, label):
        """
        Return the sequences.fasta records of gene_id.

        Each record is a dict with 'id', 'desc' and 'seq' (the fields of the
        REST answer that transcript_query writes) and 'label'. Raise KeyError
        if the gene is not in the GTF.
        """
        gene = self._con.execute(
            'SELECT seq, start, end FROM gene WHERE gene_id = ?', (gene_id,)).fetchone()
        if gene is None:
            raise KeyError(f'gene {gene_id} is not in the GTF')
        region, gene_start, gene_end = gene
        records = []
        for exon_id, start, end, strand in self._con.execute(
                'SELECT exon_id, start, end, strand FROM exon '
                'WHERE seq = ? AND start <= ? AND end >= ? ORDER BY start, exon_id',
                (region, gene_end, gene_start)):
            sequence = self.genome.fetch(region, start, end)
            if strand < 0:
                sequence = reverse_complement(sequence)
            desc = f'{self.genome.cs_assembly[region]}:{region}:{start}:{end}:{strand}'
            records.append({'label': label, 'id': exon_id, 'desc': desc,
                            'seq': sequence})
        return records
