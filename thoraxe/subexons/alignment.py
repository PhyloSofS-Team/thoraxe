"""
alignment: Module to create the subexon MSA with Clustal Omega.

This module creates a MSA of subexons using Clustal Omega.
"""

import collections
import logging
import os
import platform
import re
import shutil
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio import BiopythonWarning
from recordclass import recordclass

from thoraxe import transcript_info

ORTHO_GROUP_PATTERN = re.compile("^[0-9]+_[0-9]+")

ColPattern = collections.namedtuple('ColPattern', ['pattern', 'start', 'end'])

ColCluster = recordclass(  # pylint: disable=invalid-name
    'ColCluster', ['patterns', 'consensus', 'start', 'end'])


def subexon_connectivity(subexon_table, id_column='SubexonIDCluster'):
    """
    Return a set of connected subexon tuples.

    Each tuple has two elements, 'SubexonIDCluster' (directed) pairs.
    The first subexon appears before of the second in at least one transcript.
    """
    connected_pairs = []
    col_index = subexon_table.columns.get_loc(id_column)
    for _, transcript_df in subexon_table.groupby('TranscriptIDCluster'):
        transcript = transcript_df.sort_values(by='SubexonRank',
                                               ascending=True)
        nrows = transcript.shape[0]
        if nrows > 1:
            for row_index in range(1, nrows):
                connected_pairs.append(
                    (transcript.iloc[row_index - 1, col_index],
                     transcript.iloc[row_index, col_index]))
    return set(connected_pairs)


def _create_subexon_index(subexon_table):
    """Return a pandas' DataFrame with subexon information."""
    # NOTE : SubexonID is the same for subexons with the same sequence
    # taking phases into account. Being more specific with the subset
    # columns may cause duplicated subexons in the chimeric sequence.
    subset_columns = ['SubexonID', 'GeneID']
    unique_subexons = subexon_table.drop_duplicates(subset=subset_columns)

    unique_subexons = unique_subexons.assign(Order=[
        row.SubexonCodingStart if row.Strand == 1 else (-1 *
                                                        row.SubexonCodingEnd)
        for row in unique_subexons.itertuples()
    ])

    unique_subexons = unique_subexons.sort_values(by=['Order'])

    unique_subexons = unique_subexons.loc[:, [
        'SubexonID', 'SubexonCodingStart', 'ExonProteinSequence',
        'SubexonCodingEnd', 'SubexonRank'
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
    transcript_id_columns = ['GeneID', 'TranscriptID']
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

    for row in subexon_table.itertuples():
        subexon_matrix[row.TranscriptIndex, row.SubexonIndex] = True

    return subexon_table, subexon_matrix


def _get_sequence(subexon_info,
                  subexon_column,
                  sequence_column='ExonProteinSequence'):
    """
    Return the sequence of the subexon as a string.

    This function takes the subexon_info pandas' DataFrame as input, that has
    'SubexonIndex' as the DataFrame index.

    This deletes termination codons (*).
    """
    seq = str(subexon_info.loc[subexon_column, sequence_column])
    return seq.replace('*', '')


def create_chimeric_sequences(  # pylint: disable=too-many-locals
        subexon_table,
        subexon_matrix,
        connected_subexons,
        padding='XXXXXXXXXX'):
    """
    Create chimeric sequence for Clustal Omega.

    It returns a Dict from 'GeneID' to a tuple with the chimeric
    sequence and a Dict from 'SubexonIndex' to ...
    """
    chimerics = {}
    for gene_id, gene_df in subexon_table.groupby('GeneID'):

        # DataFrame to get a subexon information using its 'SubexonIndex'
        subexon_info_cols = [
            'SubexonIndex', 'SubexonIDCluster', 'ExonProteinSequence'
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
                    subexon_info.loc[previous_subexon, 'SubexonIDCluster'],
                    subexon_info.loc[subexon, 'SubexonIDCluster']
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
    with tempfile.NamedTemporaryFile(suffix='.fasta', delete=False,
                                     mode='w') as tmp_fasta:
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
    >>> nwin or _get_wsl_name('clustalo') is None
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


def run_aligner(chimerics, output_path='alignment.fasta', aligner='clustalo'):
    """
    Run Clustal Omega in the chimeric sequences and return the output file.

    You can pass arguments using aligner (default: 'clustalo'), e.g:
        aligner='clustalo --threads=4'

    You need Clustal Omega installed to run this function. You can install it
    from: http://www.clustal.org/omega/

    If you are using Windows 10 and you have installed Clustal Omega in
    Ubuntu using the 'Windows Subsystem for Linux', you can try with the
    following options:
        aligner='ubuntu.exe -c clustalo'
        aligner='bash.exe -c clustalo'
        aligner='wsl.exe clustalo'
    """
    if len(chimerics) == 1:
        with open(output_path, 'w') as outfile:
            _print_fasta(chimerics, outfile)
        return output_path

    input_fasta = _print_temporal_fasta(chimerics)

    command = aligner.split()

    wsl = _get_wsl_name(command[0])
    is_wsl = platform.system() == 'Windows' and wsl is not None

    if 'clustalo' in command[-1]:
        command.append('--in')

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
                    ('{} not found. Please indicate a correct aligner or '
                     'install it: http://www.clustal.org/omega/'
                     ).format(aligner))
            raise err

    return output_path


def gene2species(transcript_data):
    """Return the a dict from 'GeneID' to 'Species'."""
    return pd.Series(transcript_data.Species.values,
                     index=transcript_data['GeneID']).to_dict()


def sort_species(chimerics, gene2sp, species_list=None):
    """Sort chimerics using the gene2species output and the species_list."""
    if species_list is None:
        return chimerics

    species_order = {name: i for (i, name) in enumerate(species_list)}

    return collections.OrderedDict(
        sorted(list(chimerics.items()),
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


def _get_msa_subexon_matrix(subexon_df, chimerics, msa):
    """Return the msa as a matrix of 'SubexonIDCluster'."""
    n_seq = len(msa)
    n_col = msa.get_alignment_length()
    msa_matrix = np.empty((n_seq, n_col), dtype=object)
    msa_matrix.fill('')
    index2cluster = {
        index: cluster
        for index, cluster in zip(subexon_df['SubexonIndex'],
                                  subexon_df['SubexonIDCluster'])
    }
    for seq_index in range(0, n_seq):
        _fill_msa_matrix(msa_matrix, chimerics, msa,
                         seq_index, lambda index: index2cluster[index])

    return msa_matrix


def msa_matrices(subexon_df, chimerics, msa):
    """Return msa (residues) and subexon id cluster matrices from msa."""
    msa_matrix = np.array([list(record) for record in msa])
    subexon_matrix = _get_msa_subexon_matrix(subexon_df, chimerics, msa)
    return msa_matrix, subexon_matrix


def _compare_subexons(function, msa_matrix, subexon_matrix, min_col_number,
                      *args, **kwargs):
    """
    Compare each subexon against others in the msa using function.

    function should take an empty dictionary to update, the (query) subexon id
    and the msa matrix (numpy matrix of strings/chars) of with the sequences
    to compare.

    This function iterates each subexon of the msa, and creates a sub msa that
    only contains the columns where that subexon is present and the sequences
    that have residues in those columns.

    _compare_subexons takes 4 arguments: function, msa_matrix, subexon_matrix,
    min_col_number. Other arguments and keyword arguments are passed to
    function.

    This function returns the dictionary.
    """
    result = {}
    for subexon in np.unique(subexon_matrix):
        if subexon == '':
            continue
        # np.where returns an (rows, columns) tuple for a matrix
        start, stop = np.where(
            (subexon_matrix == subexon).any(axis=0))[0][[0, -1]]
        subexon_msa = msa_matrix[:, start:stop + 1]
        seq_indexes = np.where(np.logical_not(
            (subexon_msa == '-').all(axis=1)))
        exon_msa = subexon_msa[seq_indexes]
        n_seq, n_col = exon_msa.shape
        if n_col < min_col_number or n_seq == 1:
            continue
        function(result, subexon, exon_msa, *args, **kwargs)

    return result


def _should_keep_subexon(msa_matrix, cutoff=30.0):
    """Return True if the subexon is aligned to a similar sequence."""
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


def _add_subexon_to_delete(result, subexon, subexon_msa, cutoff=30.0):
    """Update result keys with the subexon cluster id to keep."""
    if not _should_keep_subexon(subexon_msa, cutoff=cutoff):
        result[subexon] = True


def _subexons_to_delete(msa_matrix,
                        subexon_matrix,
                        cutoff=30.0,
                        min_col_number=4):
    """Return a dict from subexon cluster id to keep to percent identity."""
    return _compare_subexons(_add_subexon_to_delete,
                             msa_matrix,
                             subexon_matrix,
                             min_col_number,
                             cutoff=cutoff)


def _delete_subexons(subexons, msa_matrix, subexon_matrix):
    """Replace an subexon by '' or '-' in subexon_matrix and msa_matrix."""
    for subexon in subexons:
        mask = subexon_matrix == subexon
        subexon_matrix[mask] = ''
        msa_matrix[mask] = '-'
    return msa_matrix, subexon_matrix


def delete_subexons(subexon_df, chimerics, msa, cutoff=30.0, min_col_number=4):
    """Return the list of 'SubexonIDCluster' to delete from 'Cluster'."""
    complete_set = set([])
    msa_matrix, subexon_matrix = msa_matrices(subexon_df, chimerics, msa)
    to_delete = _subexons_to_delete(msa_matrix,
                                    subexon_matrix,
                                    cutoff=cutoff,
                                    min_col_number=min_col_number)
    while to_delete:
        subexon_ids = to_delete.keys()
        complete_set.update(subexon_ids)
        _delete_subexons(subexon_ids, msa_matrix, subexon_matrix)
        to_delete = _subexons_to_delete(msa_matrix,
                                        subexon_matrix,
                                        cutoff=cutoff,
                                        min_col_number=min_col_number)
    return complete_set


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


def _is_orthologous_exon_group(exon):
    """
    Return True if the string has the shape number_number at the beginning.

    >>> _is_orthologous_exon_group('1_1-1_0')
    True
    >>> _is_orthologous_exon_group('2_0')
    True
    >>> _is_orthologous_exon_group('')
    False
    >>> _is_orthologous_exon_group('nan')
    False
    >>> _is_orthologous_exon_group('NaN')
    False
    """
    return ORTHO_GROUP_PATTERN.match(str(exon)) is not None


def impute_missing_orthologous_exon_group(table, column='HomologousExons'):
    """
    Replace column values that do not conform to the naming by 0_number.
    """
    number = 1
    for i in table.index:
        name = table.loc[i, column]
        if not _is_orthologous_exon_group(name):
            table.at[i, column] = '0_' + str(number)
            number += 1
    return table


def _store_homologous_subexons(subexon_df, seq, subexon, gene, exon_id):
    """
    Store homologous subexon information in subexon_df.

    exon_id is the Homologous Exons ID created by Cluster_ColCluster.
    """
    seq = seq.replace('-', '')
    length = len(seq)
    query = (subexon_df['SubexonIndex'] == subexon) & (
        subexon_df['GeneID'] == gene)
    value = subexon_df.loc[query, 'HomologousExons'].unique()[0]
    if '_' in value:
        subexon_df.loc[query, 'HomologousExons'] += '/{}'.format(exon_id)
        subexon_df.loc[query, 'HomologousExonLengths'] += '/{}'.format(length)
        subexon_df.loc[query, 'HomologousExonSequences'] += '/{}'.format(seq)
    else:
        subexon_df.loc[query, 'HomologousExons'] = exon_id
        subexon_df.loc[query, 'HomologousExonLengths'] = str(length)
        subexon_df.loc[query, 'HomologousExonSequences'] = seq


def save_homologous_subexons(subexon_df, sequences, gene_ids, colclusters,
                             output_folder):
    """
    It saves the information about homologous exons.

    It takes a list of sequences, like the one returned by msa2sequences.
    Return subexon_df with the homologous subexon information.
    For each homologous subexon saves a fasta MSA in the output_folder.
    """
    cluster = str(subexon_df['Cluster'][0])

    subexon_df = subexon_df.assign(HomologousExons=cluster,
                                   HomologousExonLengths="",
                                   HomologousExonSequences="")

    for (i, colcluster) in enumerate(colclusters):
        with open(
                os.path.join(
                    output_folder,
                    'msa_homologous_exon_{}_{}.fasta'.format(cluster, i)),
                'w') as file:
            for (j, subexon) in enumerate(colcluster.consensus):
                if not np.isnan(subexon):
                    gene = gene_ids[j]
                    seq = sequences[j][colcluster.start:colcluster.end + 1]
                    file.write('>{}\n{}\n'.format(gene, seq))
                    _store_homologous_subexons(subexon_df, seq, subexon, gene,
                                               '{}_{}'.format(cluster, i))

    return subexon_df
