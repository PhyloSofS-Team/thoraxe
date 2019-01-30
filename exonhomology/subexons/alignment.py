"""
subexon_alignment: Module to create the subexon MSA with MAFFT.

This module creates a MSA of subexons using MAFFT.
"""

import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio import BiopythonWarning
from collections import OrderedDict

# List from `orthokeep` in `EnsemblRESTTranscriptQueries.py`
# Order from NCBI Taxonomy/CommonTree, but using Human as 1
SPECIES_ORDER = {
    'homo_sapiens': 1,
    'mus_musculus': 6,
    'macaca_mulatta': 3,
    'danio_rerio': 12,
    'xenopus_tropicalis': 11,
    'caenorhabditis_elegans': 14,
    'gallus_gallus': 10,
    'rattus_norvegicus': 5,
    'bos_taurus': 7,
    'monodelphis_domestica': 4,
    'ornithorhynchus_anatinus': 9,
    'drosophila_melanogaster': 13,
    'gorilla_gorilla': 2,
    'sus_scrofa': 8
}


def subexon_connectivity(subexon_table):
    """
    Return a set of connected subexon tuples.

    Each tuple has two elements, 'Subexon ID cluster' (directed) pairs.
    The first subexon appears before of the second in at least one transcript.
    """
    connected_pairs = []
    col_index = subexon_table.columns.get_loc('Subexon ID cluster')
    for _, transcript_df in subexon_table.groupby(
            'Transcript stable ID cluster'):
        transcript = transcript_df.sort_values(
            by='Subexon rank in transcript', ascending=True)
        nrows = transcript.shape[0]
        if nrows > 1:
            for row_index in range(1, nrows):
                connected_pairs.append(
                    (transcript.iloc[row_index - 1, col_index],
                     transcript.iloc[row_index, col_index]))
    return set(connected_pairs)


def _create_subexon_index(subexon_table):
    """Return a pandas' DataFrame with subexon information."""
    # NOTE : Subexon ID is the same for subexons with the same sequence
    # taking phases into account. Being more specific with the subset
    # columns may cause duplicated subexons in the chimeric sequence.
    subset_columns = ['Subexon ID', 'Gene stable ID']
    unique_subexons = subexon_table.drop_duplicates(subset=subset_columns)

    unique_subexons = unique_subexons.assign(Order=[
        row['Subexon genomic coding start'] if row['Strand'] == 1 else (
            -1 * row['Subexon genomic coding end'])
        for _, row in unique_subexons.iterrows()
    ])

    unique_subexons = unique_subexons.sort_values(by=['Order'])

    unique_subexons = unique_subexons.loc[:, [
        'Subexon ID', 'Subexon genomic coding start', 'Exon protein sequence',
        'Subexon genomic coding end', 'Subexon rank in transcript'
    ]]

    unique_subexons = unique_subexons.assign(
        SubexonIndex=list(range(0, unique_subexons.shape[0])))

    with warnings.catch_warnings():
        # Bio/Seq.py : class Seq : __hash__ : warnings.warn
        warnings.simplefilter('ignore', BiopythonWarning)
        output = subexon_table.merge(unique_subexons)

    return output


def _create_transcript_index(subexon_table):
    """Return a pandas' DataFrame with the gene and transcript ids."""
    transcript_id_columns = ['Gene stable ID', 'Transcript stable ID']
    unique_transcripts = subexon_table.drop_duplicates(
        subset=transcript_id_columns)

    unique_transcripts = unique_transcripts.loc[:, transcript_id_columns]

    unique_transcripts = unique_transcripts.assign(
        TranscriptIndex=list(range(0, unique_transcripts.shape[0])))
    output = subexon_table.merge(unique_transcripts)

    return output


def create_subexon_matrix(subexon_table):
    """Return a binary matrix showing subexon presence in transcripts."""
    # _create_subexon_index and _create_transcript_index change
    # subexon_table with the last merge before return:
    subexon_table = _create_subexon_index(subexon_table)
    subexon_table = _create_transcript_index(subexon_table)

    subexon_table = subexon_table.sort_values(
        by=['TranscriptIndex', 'SubexonIndex'])

    n_subexons = len(subexon_table['SubexonIndex'].unique())
    n_transcripts = len(subexon_table['TranscriptIndex'].unique())

    subexon_matrix = np.zeros((n_transcripts, n_subexons), dtype=np.bool)

    for _, row in subexon_table.iterrows():
        subexon_matrix[row['TranscriptIndex'], row['SubexonIndex']] = True

    return subexon_table, subexon_matrix


def _get_sequence(subexon_info,
                  subexon_column,
                  padding,
                  sequence_column='Exon protein sequence'):
    """
    Return the sequence of the subexon as a string.

    This function takes the subexon_info pandas' DataFrame as input, that has
    'SubexonIndex' as the DataFrame index.

    This replaces termination codons (*) with the padding.
    """
    seq = str(subexon_info.loc[subexon_column, sequence_column])
    return seq.replace('*', padding)


def create_chimeric_sequences(subexon_table,
                              subexon_matrix,
                              connected_subexons,
                              padding='XXXXXXXXXX'):
    """
    Create chimeric sequence for MAFFT.

    It returns a Dict from 'Gene stable ID' to a tuple with the chimeric
    sequence and a Dict from 'SubexonIndex' to ...
    """
    chimerics = {}
    for gene_id, gene_df in subexon_table.groupby('Gene stable ID'):

        # DataFrame to get a subexon information using its 'SubexonIndex'
        subexon_info_cols = [
            'SubexonIndex', 'Subexon ID cluster', 'Exon protein sequence'
        ]
        subexon_info = gene_df.loc[:, subexon_info_cols]
        # NOTE: It make copies to not delete subexons inplace:
        with warnings.catch_warnings():
            # Bio/Seq.py : class Seq : __hash__ : warnings.warn
            warnings.simplefilter('ignore', BiopythonWarning)
            subexon_info = subexon_info.drop_duplicates(
                subset=subexon_info_cols)
        subexon_info = subexon_info.set_index('SubexonIndex')
        subexon_info = subexon_info.sort_index()

        transcript_index = sorted(gene_df['TranscriptIndex'].unique())

        subexon_index = gene_df['SubexonIndex'].unique()
        subexon_index.sort()

        transcript_matrix = subexon_matrix[transcript_index, :]
        transcript_matrix = transcript_matrix[:, subexon_index]

        subexon = subexon_index[0]
        chimeric = _get_sequence(subexon_info, subexon, padding)
        breaks = {subexon: len(chimeric)}
        for idx in range(1, len(subexon_index)):
            previous_subexon = subexon_index[idx - 1]
            subexon = subexon_index[idx]
            # Does it need padding ?
            if ((subexon_info.loc[previous_subexon, 'Subexon ID cluster'],
                 subexon_info.loc[subexon, 'Subexon ID cluster']
                 ) not in connected_subexons
                    and chimeric[-len(padding):] != padding):
                chimeric += padding
            chimeric += _get_sequence(subexon_info, subexon, padding)
            breaks.update({subexon: len(chimeric)})

        chimerics[gene_id] = (chimeric, breaks)

    return chimerics


def _print_temporal_fasta(chimerics):
    """Save chimeric sequences in a temporal fasta file and return its name."""
    with tempfile.NamedTemporaryFile(
            suffix='.fasta', delete=False, mode='w') as tmp_fasta:
        for (key, value) in chimerics.items():
            tmp_fasta.write('>%s\n%s\n' % (key, value[0]))

    return tmp_fasta.name


def run_mafft(chimerics,
              output_path='alignment.fasta',
              mafft_path='mafft',
              mafft_arguments=None):
    """
    Run MAFFT in the chimeric sequences and return the output file.

    mafft_arguments should be a list of str arguments for subprocess, e.g.:
    ['--maxiterate', '100']
    """
    input_fasta = _print_temporal_fasta(chimerics)

    command = [mafft_path, input_fasta]
    if mafft_arguments is not None:
        command.extend(mafft_arguments)

    with open(output_path, 'wb') as outfile:
        process = subprocess.Popen(command, stdout=subprocess.PIPE)
        for line in process.stdout:
            outfile.write(line)
        process.wait()

    return output_path


def gene2species(transcript_data):
    """Return the a dict from 'Gene stable ID' to 'Species'."""
    return pd.Series(
        transcript_data.Species.values,
        index=transcript_data['Gene stable ID']).to_dict()


def sort_species(chimerics, gene2sp, species_order=SPECIES_ORDER):
    """Sort chimerics using the output from gene2species and SPECIES_ORDER."""
    return OrderedDict(
        sorted(
            list(chimerics.items()),
            key=lambda x: species_order[gene2sp[x[0]]]))


def create_msa_matrix(chimerics, msa_file):  # pylint: disable=too-many-locals
    """
    Convert a msa_file from chimerics to a matrix.

    Each cell has the subexon number (Index) or nan for gaps and padding.
    """
    msa = AlignIO.read(msa_file, 'fasta')
    n_seq = len(msa[:, 0])
    assert n_seq > 1, 'There are few sequences in %s' % msa_file
    n_col = len(msa[0])
    assert n_col > 1, 'There are few columns in %s' % msa_file

    msa_matrix = np.zeros((n_seq, n_col))
    msa_matrix.fill(np.nan)

    gene_ids = []
    for seq_index in range(0, n_seq):
        record = msa[seq_index]
        gene_ids.append(record.id)
        subexon2len = chimerics[record.id][1]
        subexons = sorted(subexon2len, key=subexon2len.get)
        seq_len = 0
        subexon_index = 0
        for col_index in range(0, n_col):
            residue = record.seq[col_index]
            if residue != '-':
                seq_len += 1

            if seq_len > subexon2len[subexons[subexon_index]]:
                subexon_index += 1

            if residue in {'-', 'X'}:
                value = np.nan
            else:
                value = subexons[subexon_index]

            msa_matrix[seq_index, col_index] = value

    return gene_ids, msa_matrix


def minimal_regions(msa_matrix):
    """Return minimal homologous regions (columns) in msa_matrix."""
    non_full_nan_columns = np.where(~np.all(np.isnan(msa_matrix), axis=0))[0]
    msa_matrix = msa_matrix[:, non_full_nan_columns]
    msa_matrix[np.isnan(msa_matrix)] = -1.0
    n_seq, n_col = msa_matrix.shape

    column = msa_matrix[:, 0]
    clusters = [list(column)]
    for col_index in range(1, n_col):
        previous_column = msa_matrix[:, col_index - 1]
        column = msa_matrix[:, col_index]
        if np.any(column != previous_column) and not (
                sum(column) == -1 * n_seq
                or sum(previous_column) == -1 * n_seq):
            clusters.append(list(column))

    return clusters
