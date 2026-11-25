"""
tsl: the rows of tsl.csv from the core database dump.

transcript_query writes, for every queried gene, one tsl.csv row per transcript
returned by REST ``/overlap/id/<gene>?feature=transcript``: every current
transcript of any gene (any biotype) that overlaps the gene on the same
sequence region. The eight columns come from these REST fields:

================  ==========================  ==============================
column            REST field                  core database
================  ==========================  ==============================
Species           (queried species)
Name              external_name               xref.display_label
TranscriptID      id                          transcript.stable_id
Source            source                      transcript.source
ExperimentSource  logic_name                  analysis.logic_name
Biotype           biotype                     transcript.biotype
Flags             transcript_support_level    transcript_attrib (code TSL)
Version           version                     transcript.version
================  ==========================  ==============================

The core database dump is the only bulk file with all eight (the GTF has no
logic_name). Missing values are written as NA. The TSL value loses its "tsl"
prefix ("tsl1 (assigned to previous version 8)" becomes
"1 (assigned to previous version 8)"), as in REST. Rows are sorted by the
internal transcript_id, the order of the REST answer.
"""

import bisect
import gzip
import os
import re
from collections import defaultdict

NA = 'NA'
TABLES = ('analysis', 'attrib_type', 'transcript_attrib', 'transcript', 'xref',
          'gene')
TRANSCRIPT_COLUMNS = ['db_id', 'TranscriptID', 'Version', 'gene_stable_id',
                      'seq_region_id', 'start', 'end', 'Name', 'Source',
                      'ExperimentSource', 'Biotype', 'Flags']


def _schema_columns(sql_gz, table):
    "Column names of table, read from the CREATE TABLE of the dump schema."
    columns, inside = [], False
    with gzip.open(sql_gz, 'rt', encoding='latin1') as stream:
        for line in stream:
            if line.startswith(f'CREATE TABLE `{table}`'):
                inside = True
                continue
            if inside:
                if line.startswith(') ENGINE'):
                    break
                match = re.match(r'\s+`([^`]+)`', line)
                if match:
                    columns.append(match.group(1))
    if not columns:
        raise ValueError(f'no schema for table {table} in {sql_gz}')
    return columns


def _rows(core_dir, table, wanted=None):
    "Yield the rows of a dumped table as dicts (\\N becomes None)."
    sql_gz = os.path.join(core_dir, os.path.basename(core_dir) + '.sql.gz')
    columns = _schema_columns(sql_gz, table)
    keep = [(i, col) for i, col in enumerate(columns)
            if wanted is None or col in wanted]
    with gzip.open(os.path.join(core_dir, f'{table}.txt.gz'), 'rt',
                   encoding='latin1') as stream:
        for line in stream:
            parts = line.rstrip('\n').split('\t')
            yield {col: (None if parts[i] == r'\N' else parts[i])
                   for i, col in keep}


def _na(value):
    return NA if value is None else value


def build_tables(core_dir, out_dir):  # pylint: disable=too-many-locals
    """
    Write genes.tsv and transcripts.tsv of one species into out_dir.

    genes.tsv: gene_stable_id, seq_region_id, start, end (current genes).
    transcripts.tsv: TRANSCRIPT_COLUMNS (current transcripts), where db_id is
    the internal transcript_id.
    """
    analysis = {row['analysis_id']: row['logic_name'] for row in
                _rows(core_dir, 'analysis', {'analysis_id', 'logic_name'})}

    tsl_type = None
    for row in _rows(core_dir, 'attrib_type', {'attrib_type_id', 'code'}):
        if row['code'] == 'TSL':
            tsl_type = row['attrib_type_id']
    if tsl_type is None:
        raise ValueError(f'no TSL attrib_type in {core_dir}')
    tsl = {}
    for row in _rows(core_dir, 'transcript_attrib',
                     {'transcript_id', 'attrib_type_id', 'value'}):
        if row['attrib_type_id'] == tsl_type:
            value = row['value']
            if value.startswith('tsl'):
                value = value[len('tsl'):]
            if row['transcript_id'] in tsl:
                raise ValueError(
                    f'two TSL values for transcript {row["transcript_id"]}')
            tsl[row['transcript_id']] = value

    transcripts = [row for row in _rows(core_dir, 'transcript')
                   if row['is_current'] == '1']
    xref_ids = {row['display_xref_id'] for row in transcripts
                if row['display_xref_id']}
    xref = {row['xref_id']: row['display_label'] for row in
            _rows(core_dir, 'xref', {'xref_id', 'display_label'})
            if row['xref_id'] in xref_ids}

    genes = [row for row in _rows(core_dir, 'gene') if row['is_current'] == '1']
    gene_stable_id = {row['gene_id']: row['stable_id'] for row in genes}

    os.makedirs(out_dir, exist_ok=True)
    genes_tsv = os.path.join(out_dir, 'genes.tsv')
    with open(genes_tsv + '.part', 'w', encoding='utf-8') as stream:
        stream.write('gene_stable_id\tseq_region_id\tstart\tend\n')
        for gene in genes:
            stream.write(f"{_na(gene['stable_id'])}\t{gene['seq_region_id']}\t"
                         f"{gene['seq_region_start']}\t{gene['seq_region_end']}\n")
    transcripts_tsv = os.path.join(out_dir, 'transcripts.tsv')
    with open(transcripts_tsv + '.part', 'w', encoding='utf-8') as stream:
        stream.write('\t'.join(TRANSCRIPT_COLUMNS) + '\n')
        for row in transcripts:
            name = xref.get(row['display_xref_id']) if row['display_xref_id'] else None
            stream.write('\t'.join([
                row['transcript_id'], _na(row['stable_id']), _na(row['version']),
                _na(gene_stable_id.get(row['gene_id'])), row['seq_region_id'],
                row['seq_region_start'], row['seq_region_end'], _na(name),
                _na(row['source']), _na(analysis.get(row['analysis_id'])),
                _na(row['biotype']), tsl.get(row['transcript_id'], NA)]) + '\n')
    os.replace(genes_tsv + '.part', genes_tsv)
    # transcripts.tsv last: its presence marks the tables as complete
    os.replace(transcripts_tsv + '.part', transcripts_tsv)
    return genes_tsv, transcripts_tsv


class SpeciesTables:  # pylint: disable=too-few-public-methods
    "The tables of one species, loaded for overlap queries."

    def __init__(self, tables_dir):
        self.genes = {}
        with open(os.path.join(tables_dir, 'genes.tsv'), encoding='utf-8') as stream:
            next(stream)
            for line in stream:
                gene, region, start, end = line.rstrip('\n').split('\t')
                self.genes.setdefault(gene, (region, int(start), int(end)))
        by_region = defaultdict(list)
        with open(os.path.join(tables_dir, 'transcripts.tsv'),
                  encoding='utf-8') as stream:
            header = next(stream).rstrip('\n').split('\t')
            for line in stream:
                record = dict(zip(header, line.rstrip('\n').split('\t')))
                record['start'], record['end'] = int(record['start']), int(record['end'])
                by_region[record['seq_region_id']].append(record)
        self.regions = {}
        for region, records in by_region.items():
            records.sort(key=lambda record: record['start'])
            starts = [record['start'] for record in records]
            longest = max(record['end'] - record['start'] for record in records)
            self.regions[region] = (records, starts, longest)

    def rows(self, species, gene_id):
        """
        Return the tsl.csv rows of gene_id as dicts (REST field names).

        Every current transcript overlapping the gene on its sequence region,
        in internal transcript_id order. Empty if gene_id is not a current gene.
        """
        if gene_id not in self.genes:
            return []
        region, gene_start, gene_end = self.genes[gene_id]
        if region not in self.regions:
            return []
        records, starts, longest = self.regions[region]
        first = bisect.bisect_left(starts, gene_start - longest)
        last = bisect.bisect_right(starts, gene_end)
        hits = [record for record in records[first:last]
                if record['end'] >= gene_start]
        hits.sort(key=lambda record: int(record['db_id']))
        return [{'species': species,
                 'external_name': record['Name'],
                 'transcript_id': record['TranscriptID'],
                 'source': record['Source'],
                 'logic_name': record['ExperimentSource'],
                 'biotype': record['Biotype'],
                 'transcript_support_level': record['Flags'],
                 'version': record['Version']} for record in hits]


class TablesCache:
    """
    Keep the tables of the last `size` species in memory.

    The query species and one ortholog species at a time is enough: the
    orthologs of a gene are queried species by species.
    """

    def __init__(self, bulk, size=2):
        self.bulk = bulk
        self.size = size
        self._cache = {}

    def get(self, species):
        "Return the SpeciesTables of species."
        if species in self._cache:
            self._cache[species] = self._cache.pop(species)
        else:
            self._cache[species] = SpeciesTables(self.bulk.tsl_tables_dir(species))
            while len(self._cache) > self.size:
                self._cache.pop(next(iter(self._cache)))
        return self._cache[species]

    def rows(self, species, gene_id):
        "Return the tsl.csv rows of gene_id in species."
        return self.get(species).rows(species, gene_id)
