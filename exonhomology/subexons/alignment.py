"""
alignment: Module to create the subexon MSA with MAFFT.

This module creates a MSA of subexons using MAFFT.
"""

import collections
import logging
import os
import platform
import shutil
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio import BiopythonWarning
from recordclass import recordclass

from exonhomology import transcript_info

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

ColPattern = collections.namedtuple('ColPattern', ['pattern', 'start', 'end'])

ColCluster = recordclass(  # pylint: disable=invalid-name
    'ColCluster', ['patterns', 'consensus', 'start', 'end'])


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
                  sequence_column='Exon protein sequence'):
    """
    Return the sequence of the subexon as a string.

    This function takes the subexon_info pandas' DataFrame as input, that has
    'SubexonIndex' as the DataFrame index.

    This deletes termination codons (*).
    """
    seq = str(subexon_info.loc[subexon_column, sequence_column])
    return seq.replace('*', '')


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
        chimeric = _get_sequence(subexon_info, subexon)
        breaks = {subexon: len(chimeric)}
        for idx in range(1, len(subexon_index)):
            previous_subexon = subexon_index[idx - 1]
            subexon = subexon_index[idx]
            # Does it need padding ?
            if chimeric and not chimeric.endswith(padding) and (
                    subexon_info.loc[previous_subexon, 'Subexon ID cluster'],
                    subexon_info.loc[subexon, 'Subexon ID cluster']
            ) not in connected_subexons:
                chimeric += padding
            chimeric += _get_sequence(subexon_info, subexon)
            breaks.update({subexon: len(chimeric)})

        if chimeric:
            chimerics[gene_id] = (chimeric, breaks)

    return chimerics


def _print_fasta(chimerics, stream):
    """Write chimeric sequences in fasta format."""
    for (key, value) in chimerics.items():
        chimeric = value[0]
        if chimeric:
            stream.write('>{}\n{}\n'.format(key, chimeric))


def _print_temporal_fasta(chimerics):
    """Save chimeric sequences in a temporal fasta file and return its name."""
    with tempfile.NamedTemporaryFile(
            suffix='.fasta', delete=False, mode='w') as tmp_fasta:
        _print_fasta(chimerics, tmp_fasta)

    return tmp_fasta.name


def _get_wsl_name(executable_path):
    r"""
    Return the name of the 'Windows Subsystem for Linux' executable.

    Return None is the path in not valid.

    >>> import platform
    >>> nwin = platform.system() != 'Windows'
    >>> nwin or _get_wsl_name('bash.exe') == 'bash.exe'
    True
    >>> nwin or _get_wsl_name('C:\\WINDOWS\\SysNative\\bash.exe') == 'bash.exe'
    True
    >>> nwin or _get_wsl_name('mafft') is None
    True
    """
    executable_name = os.path.basename(os.path.abspath(executable_path))
    return executable_name if executable_name.split('.')[0].lower() in {
        'wsl', 'ubuntu', 'bash'
    } else None


def _get_wsl_path(executable_name):
    """
    Return path to the 'Windows Subsystem for Linux' executable.

    Stackoverflow: python-subprocess-call-cannot-find-windows-bash-exe
    """
    is32bit = platform.architecture()[0] == '32bit'
    system32 = os.path.join(os.environ['SystemRoot'],
                            'SysNative' if is32bit else 'System32')
    return os.path.join(system32, executable_name)


def _win2wsl(path):
    r"""
    Convert a Windows path to a 'Windows Subsystem for Linux' path.

    This is similar to wslpath.
    Stackoverflow: python-subprocess-call-cannot-find-windows-bash-exe

    >>> import platform
    >>> nowindows = platform.system() != 'Windows'
    >>> nowindows or _win2wsl('C:\\aa\\bb\\foo.zip') == '/mnt/c/aa/bb/foo.zip'
    True
    """
    path = os.path.abspath(path)
    if len(path) > 3 and path[1:2] == ':':
        drive = path[:1].lower()
        return '/mnt/' + drive + path[2:].replace('\\', '/')

    raise ValueError(
        '{} is not an absolute Windows path to a file.'.format(path))


def run_mafft(chimerics, output_path='alignment.fasta', mafft_path='mafft'):
    """
    Run MAFFT in the chimeric sequences and return the output file.

    You can pass arguments using mafft_path (default: 'mafft'), e.g:
        mafft_path='mafft --maxiterate 100 --auto'

    You need MAFFT installed to run this function. You can install MAFFT from:
        https://mafft.cbrc.jp/alignment/software/

    If you are using Windows 10 and you have installed MAFFT in Ubuntu using
    the 'Windows Subsystem for Linux', you can try with the following options:
        mafft_path='ubuntu.exe -c mafft'
        mafft_path='bash.exe -c mafft'
        mafft_path='wsl.exe mafft'
    """
    if len(chimerics) == 1:
        with open(output_path, 'w') as outfile:
            _print_fasta(chimerics, outfile)
        return output_path

    input_fasta = _print_temporal_fasta(chimerics)

    command = mafft_path.split()

    wsl = _get_wsl_name(command[0])
    is_wsl = platform.system() == 'Windows' and wsl is not None

    if is_wsl:
        command[0] = _get_wsl_path(wsl)
        command.append(_win2wsl(input_fasta))
        command.append('>')
        command.append(_win2wsl(output_path))
        if wsl.lower().startswith('wsl'):
            subprocess.call(command)
        else:
            assert command[1] == '-c'
            subprocess.call("{} -c '{}'".format(command[0],
                                                ' '.join(command[2:])))
    else:
        try:
            command.append(input_fasta)
            with open(output_path, 'wb') as outfile:
                process = subprocess.Popen(command, stdout=subprocess.PIPE)
                for line in process.stdout:
                    outfile.write(line)
                process.wait()
        except (OSError, FileNotFoundError) as err:
            if shutil.which(command[0]) is None:
                raise OSError(
                    ('{} not found. Please indicate a correct mafft_path or '
                     'install it: https://mafft.cbrc.jp/alignment/software/'
                     ).format(mafft_path))
            else:
                raise err

    return output_path


def gene2species(transcript_data):
    """Return the a dict from 'Gene stable ID' to 'Species'."""
    return pd.Series(
        transcript_data.Species.values,
        index=transcript_data['Gene stable ID']).to_dict()


def sort_species(chimerics, gene2sp, species_order=None):
    """Sort chimerics using the output from gene2species and SPECIES_ORDER."""
    if species_order is None:
        species_order = SPECIES_ORDER

    return collections.OrderedDict(
        sorted(
            list(chimerics.items()),
            key=lambda x: species_order[gene2sp[x[0]]]))


def read_msa_fasta(msa_file):
    """
    Return a BioPython's alignment object from the fasta msa_file.

    Return None if the file is empty.
    """
    if os.path.getsize(msa_file) > 0:
        msa = AlignIO.read(msa_file, 'fasta')
        n_seq = len(msa)
        if n_seq < 2:
            logging.warning('There are few (%s) sequences in %s', n_seq,
                            msa_file)
        n_col = msa.get_alignment_length()
        if n_col < 2:
            logging.warning('There are few (%s) columns in %s', n_col,
                            msa_file)
        return msa
    return None


def get_gene_ids(msa):
    """Take a biopython MSA and return a list with the gene (sequence) ids."""
    return [record.id for record in msa]


def _fill_msa_matrix(msa_matrix,
                     chimerics,
                     msa,
                     seq_index,
                     transform_subexon_index=lambda index: index):
    """Function to fill msa_matrix with msa data."""
    record = msa[seq_index]
    subexon2len = chimerics[record.id][1]
    subexons = sorted(subexon2len, key=subexon2len.get)
    seq_len = 0
    subexon_index = 0
    for col_index in range(0, msa_matrix.shape[1]):
        residue = record.seq[col_index]

        if residue != '-':
            seq_len += 1

        if seq_len > subexon2len[subexons[subexon_index]]:
            subexon_index += 1

        if residue not in {'-', 'X'}:
            msa_matrix[seq_index, col_index] = transform_subexon_index(
                subexons[subexon_index])


def _get_msa_exon_matrix(subexon_df, chimerics, msa):
    """Return the msa as a matrix of 'Exon stable ID cluster'."""
    n_seq = len(msa)
    n_col = msa.get_alignment_length()
    msa_matrix = np.empty((n_seq, n_col), dtype=object)
    msa_matrix.fill('')
    index2cluster = {
        index: cluster
        for index, cluster in zip(subexon_df['SubexonIndex'],
                                  subexon_df['Exon stable ID cluster'])
    }
    for seq_index in range(0, n_seq):
        _fill_msa_matrix(msa_matrix, chimerics, msa,
                         seq_index, lambda index: index2cluster[index])

    return msa_matrix


def msa_matrices(subexon_df, chimerics, msa):
    """Return msa (residues) and exon (exon cluster id) matrices from msa."""
    msa_matrix = np.array([list(record) for record in msa])
    exon_matrix = _get_msa_exon_matrix(subexon_df, chimerics, msa)
    return msa_matrix, exon_matrix


def compare_exons(msa_matrix, exon_matrix, function, *args, **kwargs):
    """
    Compare each exon against others in the msa using function.

    function should take an empty dictionary to update, the (query) exon id
    and the msa matrix (numpy matrix of strings/chars) of with the sequences
    to compare.

    This function iterates each exon of the msa, and creates a sub msa that
    only contains the columns where that exon is present and the sequences
    that have residues in those columns.

    compare_exons takes 3 arguments: msa_matrix, exon_matrix, function.
    Other arguments and keyword arguments are passed to function.

    This function returns the dictionary.
    """
    result = {}
    for exon in np.unique(exon_matrix):
        if exon == '':
            continue
        # np.where returns an (rows, columns) tuple for a matrix
        start, stop = np.where((exon_matrix == exon).any(axis=0))[0][[0, -1]]
        exon_msa = msa_matrix[:, start:stop + 1]
        seq_indexes = np.where(np.logical_not((exon_msa == '-').all(axis=1)))
        exon_msa = exon_msa[seq_indexes]
        function(result, exon, exon_msa, *args, **kwargs)

    return result


def _should_keep_exon(msa_matrix, cutoff=30.0):
    n_seqs = msa_matrix.shape[0]
    if n_seqs == 1:
        return False
    for i in range(0, n_seqs - 1):
        query = msa_matrix[i, :]
        for j in range(i + 1, n_seqs):
            target = msa_matrix[j, :]
            pid = transcript_info.percent_identity(query, target)
            if pid >= cutoff:
                return True
    return False


def _add_exon_to_delete(result, exon, exon_msa, cutoff=30.0):
    """Update result keys with the exon cluster id to keep."""
    if not _should_keep_exon(exon_msa, cutoff=cutoff):
        result[exon] = True


def exons_to_delete(msa_matrix, exon_matrix, cutoff=30.0):
    """Return a dict from exon cluster id to keep to percent identity."""
    return compare_exons(
        msa_matrix, exon_matrix, _add_exon_to_delete, cutoff=cutoff)


def delete_exons(exons, msa_matrix, exon_matrix):
    """Replace an exon by '' or '-' in exon_matrix and msa_matrix."""
    for exon in exons:
        mask = exon_matrix == exon
        exon_matrix[mask] = ''
        msa_matrix[mask] = '-'
    return msa_matrix, exon_matrix


def create_msa_matrix(chimerics, msa):
    """
    Convert a msa from chimerics to a matrix.

    Each cell has the subexon number (Index) or nan for gaps and padding.
    """
    n_seq = len(msa)
    n_col = msa.get_alignment_length()
    msa_matrix = np.zeros((n_seq, n_col))
    msa_matrix.fill(np.nan)
    for seq_index in range(0, n_seq):
        _fill_msa_matrix(msa_matrix, chimerics, msa, seq_index)

    return msa_matrix


def column_patterns(msa_matrix):
    """Return a ColPattern list from the msa_matrix."""
    colpatterns = []
    for i in range(msa_matrix.shape[1]):
        col = msa_matrix[:, i].copy()
        if i == 0:
            colpatterns.append(ColPattern(col, i, i))
        else:
            previous_pattern = colpatterns[-1]
            if np.allclose(col, previous_pattern.pattern, equal_nan=True):
                colpatterns[-1] = ColPattern(col, previous_pattern.start, i)
            else:
                colpatterns.append(ColPattern(col, i, i))
    return colpatterns


def _equal_without_nans(col_i, col_j):
    """
    Return True if the columns are equal without comparing rows with nans.

    Returns False otherwise or if there are no elements in common between the
    columns.
    """
    nans = np.isnan(col_i) | np.isnan(col_j)
    elements_i = col_i[~nans]
    elements_j = col_j[~nans]
    if elements_i.size > 0:
        return np.array_equal(elements_i, elements_j)
    return False


def _colcluster(colpattern):
    """Return a ColCluster with colpattern as the unique element."""
    return ColCluster([colpattern], colpattern.pattern.copy(),
                      colpattern.start, colpattern.end)


def column_clusters(colpatterns):
    """Return a ColCluster liste from a ColPattern list."""
    n_patterns = len(colpatterns)
    colpattern = colpatterns[0]
    clusters = [_colcluster(colpattern)]
    if n_patterns > 1:
        for i in range(1, n_patterns):
            colpattern = colpatterns[i]
            if not np.all(np.isnan(colpattern.pattern)):
                cluster = clusters[-1]
                if (colpattern.start - cluster.end == 1
                        and _equal_without_nans(cluster.consensus,
                                                colpattern.pattern)):
                    nans = np.isnan(cluster.consensus)
                    cluster.consensus[nans] = colpattern.pattern[nans]
                    cluster.patterns.append(colpattern)
                    cluster.end = colpattern.end
                else:
                    clusters.append(_colcluster(colpattern))
    return clusters


def msa2sequences(msa, gene_ids, padding):
    """
    Return str sequences from msa.

    It also checks gene_ids and replaces padding by gaps.
    """
    sequences = []
    for i, seq in enumerate(msa):
        assert seq.id == gene_ids[i]
        sequences.append(str(seq.seq).replace(padding, '-' * len(padding)))
    return sequences


def _store_homologous_subexons(subexon_df, seq, subexon, gene, exon_id):
    """
    Store homologous exon information in subexon_df.

    exon_id is the Homologous Exons ID created by Cluster_ColCluster.
    """
    seq = seq.replace('-', '')
    length = len(seq)
    query = (subexon_df['SubexonIndex'] == subexon) & (
        subexon_df['Gene stable ID'] == gene)
    value = subexon_df.loc[query, 'HomologousExons'].unique()[0]
    if '_' in value:
        subexon_df.loc[query, 'HomologousExons'] += '-{}'.format(exon_id)
        subexon_df.loc[query, 'HomologousExonLengths'] += '-{}'.format(length)
        subexon_df.loc[query, 'HomologousExonSequences'] += '-{}'.format(seq)
    else:
        subexon_df.loc[query, 'HomologousExons'] = exon_id
        subexon_df.loc[query, 'HomologousExonLengths'] = str(length)
        subexon_df.loc[query, 'HomologousExonSequences'] = seq


def save_homologous_subexons(subexon_df, sequences, gene_ids, colclusters,
                             output_folder):
    """
    It saves the information about homologous exons.

    It takes a list of sequences, like the one returned by msa2sequences.
    Return subexon_df with the homologous exon information.
    For each homologous exon saves a fasta MSA in the output_folder.
    """
    cluster = str(subexon_df['Cluster'][0])
    subexon_df = subexon_df.assign(
        HomologousExons=cluster,
        HomologousExonLengths="",
        HomologousExonSequences="")
    for (i, colcluster) in enumerate(colclusters):
        with open(
                os.path.join(
                    output_folder, 'msa_homologous_exon_{}_{}.fasta'.format(
                        cluster, i)), 'w') as file:
            for (j, subexon) in enumerate(colcluster.consensus):
                if not np.isnan(subexon):
                    gene = gene_ids[j]
                    seq = sequences[j][colcluster.start:colcluster.end + 1]
                    file.write('>{}\n{}\n'.format(gene, seq))
                    _store_homologous_subexons(subexon_df, seq, subexon, gene,
                                               '{}_{}'.format(cluster, i))

    return subexon_df
