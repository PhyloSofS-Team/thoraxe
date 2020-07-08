"""
alignment: Module to create the subexon MSA with ProGraphMSA.
"""

# pylint: disable=too-many-lines

import collections
import logging
import os
import platform
import re
import shutil
import subprocess
from tempfile import NamedTemporaryFile
import warnings

import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio import BiopythonWarning
from recordclass import recordclass

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
                    (transcript.iloc[row_index - 1,
                                     col_index], transcript.iloc[row_index,
                                                                 col_index]))
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
        'SubexonID', 'SubexonCodingStart', 'SubexonProteinSequence',
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
                  sequence_column='SubexonProteinSequence'):
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
    Create chimeric sequence for ProGraphMSA.

    It returns a Dict from 'GeneID' to a tuple with the chimeric
    sequence and a Dict from 'SubexonIndex' to ...
    """
    chimerics = {}
    for gene_id, gene_df in subexon_table.groupby('GeneID'):

        # DataFrame to get a subexon information using its 'SubexonIndex'
        subexon_info_cols = [
            'SubexonIndex', 'SubexonIDCluster', 'SubexonProteinSequence'
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
    with NamedTemporaryFile(suffix='.fasta', delete=False, mode='w') as tmp:
        _print_fasta(chimerics, tmp)

    return tmp.name


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
    >>> nwin or _get_wsl_name('ProGraphMSA') is None
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


def run_aligner(chimerics,
                output_path='alignment.fasta',
                aligner='ProGraphMSA'):
    """
    Run ProGraphMSA in the chimeric sequences and return the output file.

    You can pass arguments using aligner (default: 'ProGraphMSA'), e.g:
    ::

        aligner='ProGraphMSA --mldist_gap'

    You need ProGraphMSA installed to run this function. You can install it
    from: https://github.com/acg-team/ProGraphMSA

    If you are using Windows 10 and you have installed ProGraphMSA in
    Ubuntu using the 'Windows Subsystem for Linux', you can try with the
    following options:
    ::

        aligner='ubuntu.exe -c ProGraphMSA'
        aligner='bash.exe -c ProGraphMSA'
        aligner='wsl.exe ProGraphMSA'
    """
    if len(chimerics) == 1:
        with open(output_path, 'w') as outfile:
            _print_fasta(chimerics, outfile)
        return output_path

    input_fasta = _print_temporal_fasta(chimerics)

    command = aligner.split()

    wsl = _get_wsl_name(command[0])
    is_wsl = platform.system() == 'Windows' and wsl is not None

    if '--input_order' not in command:
        command.append('--input_order')
    if '--fasta' not in command:
        command.append('--fasta')

    if is_wsl:
        command[0] = _get_wsl_path(wsl)
        command.append('--output')
        command.append(_win2wsl(output_path))
        command.append(_win2wsl(input_fasta))
        if wsl.lower().startswith('wsl'):
            subprocess.call(command)
        else:
            assert command[1] == '-c'
            subprocess.call("{} -c '{}'".format(command[0],
                                                ' '.join(command[2:])))
    else:
        try:
            command.append('--output')
            command.append(output_path)
            command.append(input_fasta)
            subprocess.call(command)
        except (OSError, FileNotFoundError) as err:
            if shutil.which(command[0]) is None:
                raise OSError(
                    '{} not found. '.format(aligner) +
                    'Please indicate the path to ProGraphMSA or install it ' +
                    'from https://github.com/acg-team/ProGraphMSA')
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
    if not chimerics:
        return np.empty((0, 0), dtype=object)
    n_seq = len(msa)
    n_col = msa.get_alignment_length()
    msa_matrix = np.empty((n_seq, n_col), dtype=object)
    msa_matrix.fill('')
    index2cluster = dict(
        zip(subexon_df['SubexonIndex'], subexon_df['SubexonIDCluster']))
    for seq_index in range(0, n_seq):
        _fill_msa_matrix(msa_matrix, chimerics, msa, seq_index,
                         lambda index: index2cluster[index])

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


def _percent_identity_without_gaps(query, target):
    """
    Percent identity of aligned sequences without counting gaps as mismatches.

    >>> _percent_identity_without_gaps("AA---", "GAAA-")
    50.0
    """
    aln_len = len(query)
    identical = 0
    for res_i, res_j in zip(query, target):
        if (res_i == '-') or (res_j == '-'):
            aln_len -= 1
        elif res_i == res_j:
            identical += 1
    if not aln_len:
        return np.nan
    return 100.0 * (identical / aln_len)


def _should_keep_subexon(msa_matrix, cutoff=30.0, keep_single_subexons=False):
    """Return True if the subexon is aligned to a similar sequence."""
    n_seqs = msa_matrix.shape[0]
    if n_seqs == 1:
        return keep_single_subexons
    tries = 0
    for i in range(0, n_seqs - 1):
        query = msa_matrix[i, :]
        for j in range(i + 1, n_seqs):
            target = msa_matrix[j, :]
            pid = _percent_identity_without_gaps(query, target)
            if np.isnan(pid):
                continue
            if pid >= cutoff:
                return True
            tries += 1
    if not tries:
        # keep_single_subexons=True : do not migrate/discard sub-exons that
        # align only with gaps.
        return keep_single_subexons

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
    if not chimerics:
        return np.empty((0, 0), dtype=object)
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
    columns. Also, it returns True if col_i or col_j are only gap/padding
    columns, i.e. full nan columns.

    >>> import numpy as np
    >>> _equal_without_nans(np.array([1, 2, 3]), np.array([1, 2, 3]))
    True
    >>> _equal_without_nans(np.array([1, np.nan, 3]), np.array([1, 2, np.nan]))
    True
    >>> _equal_without_nans(np.array([1, 2]), np.array([np.nan, np.nan]))
    True
    >>> _equal_without_nans(np.array([np.nan, np.nan]), np.array([1, 2]))
    True
    >>> _equal_without_nans(np.array([0, 2, 3]), np.array([1, 2, 3]))
    False
    >>> _equal_without_nans(np.array([0, np.nan, 3]), np.array([1, 2, np.nan]))
    False
    """
    nan_i = np.isnan(col_i)
    nan_j = np.isnan(col_j)
    if np.all(nan_i) or np.all(nan_j):
        return True
    nan_rows = nan_i | nan_j
    elements_i = col_i[~nan_rows]
    elements_j = col_j[~nan_rows]
    if elements_i.size > 0:
        return np.array_equal(elements_i, elements_j)
    return False


def _colcluster(colpattern):
    """Return a ColCluster with colpattern as the unique element."""
    return ColCluster([colpattern], colpattern.pattern.copy(),
                      colpattern.start, colpattern.end)


def _cluster_column_clusters(col_clusters):
    """
    Cluster column clusters that do not change sub-exons in their consensus.
    """
    n_clusters = len(col_clusters)
    if not n_clusters:
        return []
    merged_clusters = [col_clusters[0]]
    if n_clusters > 1:
        for i in range(1, n_clusters):
            current_cluster = col_clusters[i]
            last_cluster = merged_clusters[-1]
            if (current_cluster.start - last_cluster.end == 1
                    and _equal_without_nans(last_cluster.consensus,
                                            current_cluster.consensus)):
                nans = np.isnan(last_cluster.consensus)
                last_cluster.consensus[nans] = current_cluster.consensus[nans]
                last_cluster.patterns.extend(current_cluster.patterns)
                last_cluster.end = current_cluster.end
            else:
                merged_clusters.append(current_cluster)
    return merged_clusters


def column_clusters(colpatterns):
    """Return a ColCluster list from a ColPattern list."""
    n_patterns = len(colpatterns)
    if not n_patterns:
        return []
    colpattern = colpatterns[0]
    clusters = [_colcluster(colpattern)]
    if n_patterns > 1:
        for i in range(1, n_patterns):
            colpattern = colpatterns[i]
            cluster = clusters[-1]
            if (colpattern.start - cluster.end == 1 and _equal_without_nans(
                    cluster.consensus, colpattern.pattern)):
                nans = np.isnan(cluster.consensus)
                cluster.consensus[nans] = colpattern.pattern[nans]
                cluster.patterns.append(colpattern)
                cluster.end = colpattern.end
            else:
                clusters.append(_colcluster(colpattern))
    return _cluster_column_clusters(clusters)


def _delete_padding(seq, padding):
    """
    Replace padding by gaps.

    >>> _delete_padding("MHGL--XXXX-XXX--XXXKLMHSXXXXXXX-X-X-X", "XXXXXXXXXX")
    'MHGL---------------KLMHS-------------'

    >>> _delete_padding("MHGLXKLMHSX", "XXXXXXXXXX")
    'MHGLXKLMHSX'
    """
    regex = ''.join(c + '-*' for c in padding)
    match = re.search(regex, seq)
    if match is not None:
        actual_padding = match.group()
        seq = seq.replace(actual_padding, '-' * len(actual_padding))
        return _delete_padding(seq, padding)

    return seq


def msa2sequences(msa, gene_ids, padding):
    """
    Return str sequences from msa.

    It also checks gene_ids and replaces padding by gaps.
    Full gap columns are deleted.
    """
    sequences = []
    for i, seq in enumerate(msa):
        assert seq.id == gene_ids[i]
        sequences.append(_delete_padding(str(seq.seq), padding))
    return sequences


def _is_s_exon(exon):
    """
    Return True if the string has the shape number_number at the beginning.

    >>> _is_s_exon('1_1-1_0')
    True
    >>> _is_s_exon('2_0')
    True
    >>> _is_s_exon('')
    False
    >>> _is_s_exon('nan')
    False
    >>> _is_s_exon('NaN')
    False
    """
    return ORTHO_GROUP_PATTERN.match(str(exon)) is not None


def impute_missing_s_exon(table, column='S_exons'):
    """
    Replace column values that do not conform to the naming by 0_number.
    """
    number = 1
    for i in table.index:
        name = table.loc[i, column]
        if not _is_s_exon(name):
            table.at[i, column] = '0_' + str(number)
            number += 1
    return table


def _store_s_exons(subexon_df, seq, subexon, gene, s_exon_id):
    """
    Store s-exon information in subexon_df.

    s_exon_id is the ID of the orthologous exonic region (s-exon) created as
    Cluster_ColCluster.
    """
    seq = seq.replace('-', '')
    length = len(seq)
    query = (subexon_df['SubexonIndex'] == subexon) & (subexon_df['GeneID']
                                                       == gene)
    value = subexon_df.loc[query, 'S_exons'].unique()[0]
    if '_' in value:
        subexon_df.loc[query, 'S_exons'] += '/{}'.format(s_exon_id)
        subexon_df.loc[query, 'S_exon_Lengths'] += '/{}'.format(length)
        subexon_df.loc[query, 'S_exon_Sequences'] += '/{}'.format(seq)
    else:
        subexon_df.loc[query, 'S_exons'] = s_exon_id
        subexon_df.loc[query, 'S_exon_Lengths'] = str(length)
        subexon_df.loc[query, 'S_exon_Sequences'] = seq


def _full_gap_columns(sequences):
    """
    Return a boolean list, where True indicates a full gap column.
    """
    gaps = [residue == '-' for residue in sequences[0]]
    n_seqs = len(sequences)
    if n_seqs > 1:
        for i in range(1, n_seqs):
            sequence = sequences[i]
            for (j, residue) in enumerate(sequence):
                if residue != '-':
                    gaps[j] = False
    return gaps


def _delete_full_gap_columns(sequences):
    """
    Return the sequence without the full gap columns.
    """
    full_gaps = _full_gap_columns(sequences)
    cleaned_seqs = []
    for sequence in sequences:
        cleaned_seq = []
        for (full_gap, residue) in zip(full_gaps, sequence):
            if not full_gap:
                cleaned_seq.append(residue)
        cleaned_seqs.append(''.join(cleaned_seq))
    return cleaned_seqs


def save_s_exons(subexon_df, sequences, gene_ids, colclusters, output_folder):
    """
    It saves the information about s-exons.

    It takes a list of sequences, like the one returned by msa2sequences.
    Return subexon_df with the orthologous exonic region (s-exon) information.
    For each s-exon saves a fasta MSA in the output_folder.
    """
    if not subexon_df.empty:
        cluster = str(subexon_df['Cluster'][0])

        subexon_df = subexon_df.assign(S_exons=cluster,
                                       S_exon_Lengths="",
                                       S_exon_Sequences="")

        sequences = _delete_full_gap_columns(sequences)

        for (i, colcluster) in enumerate(colclusters):
            with open(
                    os.path.join(output_folder,
                                 'msa_s_exon_{}_{}.fasta'.format(cluster, i)),
                    'w') as file:
                for (j, subexon) in enumerate(colcluster.consensus):
                    if not np.isnan(subexon):
                        gene = gene_ids[j]
                        seq = sequences[j][colcluster.start:colcluster.end + 1]
                        file.write('>{}\n{}\n'.format(gene, seq))
                        _store_s_exons(subexon_df, seq, subexon, gene,
                                       '{}_{}'.format(cluster, i))

    return subexon_df


# Fix sub-exon blocks in the chimeric MSA


def _same_subexon(id_a, id_b):
    """
    Return True if the sub-exon ids are identical or missing (NaNs).

    >>> import numpy as np
    >>> _same_subexon(np.nan, np.nan)
    True
    >>> _same_subexon(np.nan, 2.0)
    False
    >>> _same_subexon(2.0, 2.0)
    True
    >>> _same_subexon(2.0, 8.0)
    False
    """
    if np.isnan(id_a) and np.isnan(id_b):
        return True
    return id_a == id_b


def resume_seq(matrix_row):
    """
    Return the sub-exon and gap blocks in the msa matrix.

    The function return three lists, a list of sub-exon ids and NaNs, i.e.
    the block id, a list with the column indices where the block starts and a
    list with the block ends.

    >>> resume_seq([1.0, 1.0, 1.0, np.nan, np.nan, np.nan, 1.0, 2.0, 2.0, 2.0])
    ([1.0, nan, 1.0, 2.0], [0, 3, 6, 7], [3, 6, 7, 10])
    """
    subexons = []
    starts = []
    ends = []
    for (i, subexon) in enumerate(matrix_row):
        if i == 0:
            subexons.append(subexon)
            starts.append(i)
        else:
            if not _same_subexon(subexon, subexons[-1]):
                subexons.append(subexon)
                ends.append(i)
                starts.append(i)
    ends.append(len(matrix_row))
    return (subexons, starts, ends)


ProblematicSubexonBlock = collections.namedtuple('ProblematicSubexonBlock', [
    'sequence_index', 'subexon', 'subexon_blocks', 'block_type',
    'subexon_block_start', 'subexon_block_end', 'gap_block_start',
    'gap_block_end'
])


def _is_problematic(  # pylint: disable=too-many-arguments
        starts, ends, block_index, first, max_res_block, min_gap_block):
    """
    Return True is the sub-exon block is problematic.

    In particular, a sub-exon block is problematic if it has <= max_res_block
    and it is separeated from the rest of the sub-exon by a gap block
    => min_gap_block.

    >>> _is_problematic([0, 3, 6, 7], [3, 6, 7, 10], 0, True, 2, 2)
    False
    >>> _is_problematic([0, 3, 6, 7], [3, 6, 7, 10], 2, False, 2, 2)
    True
    """
    subexon_len = ends[block_index] - starts[block_index]
    if first:
        gaps_len = ends[block_index + 1] - starts[block_index + 1]
    else:
        gaps_len = ends[block_index - 1] - starts[block_index - 1]
    return (subexon_len <= max_res_block) and (gaps_len >= min_gap_block)


def _store_problematic(  # pylint: disable=too-many-arguments
        problematic_list, sequence, subexon, subexon_blocks, starts, ends,
        first_block, last_block, max_res_block, min_gap_block):
    """
    Store problematic sub-exon blocks in the problematic_list.
    """
    if subexon_blocks > 1:
        if _is_problematic(starts, ends, first_block, True, max_res_block,
                           min_gap_block):
            problematic_list.append(
                ProblematicSubexonBlock(sequence, subexon, subexon_blocks,
                                        "first_block", starts[first_block],
                                        ends[first_block],
                                        starts[first_block + 1],
                                        ends[first_block + 1]))
        if _is_problematic(starts, ends, last_block, False, max_res_block,
                           min_gap_block):
            problematic_list.append(
                ProblematicSubexonBlock(sequence, subexon, subexon_blocks,
                                        "last_block", starts[last_block],
                                        ends[last_block],
                                        starts[last_block - 1],
                                        ends[last_block - 1]))


def _choose_one(blocks):
    """
    Search for two-blocks sub-exons that has both blocks and choose one.

    It chooses the smaller one or the first one.
    """
    subexon2blocks = collections.Counter(
        [block.subexon for block in blocks if block.subexon_blocks == 2])
    problematic = [
        subexon for (subexon, blocks) in subexon2blocks.items() if blocks == 2
    ]
    pos2len = {}
    for subexon in problematic:
        pos2len = {
            pos: block.subexon_block_end - block.subexon_block_start
            for (pos, block) in enumerate(blocks) if block.subexon == subexon
        }
        pos_to_delete = max(pos2len, key=lambda key: pos2len[key])
        del blocks[pos_to_delete]
    return blocks


def problematic_subexon_blocks(  # pylint: disable=too-many-arguments
        subexons,
        starts,
        ends,
        sequence,
        max_res_block=2,  # to move
        min_gap_block=2  # to consider
):
    """
    Return a list of ProblematicSubexonBlocks.

    In particular, a sub-exon block is problematic if it has <= max_res_block
    and it is separeated from the rest of the sub-exon by a gap block
    => min_gap_block.

    >>> import numpy as np
    >>> problematic_subexon_blocks([1.0, np.nan, 1.0, 2.0], [0, 3, 9, 10], \
[3, 9, 10, 13], 0)
    [ProblematicSubexonBlock(sequence_index=0, subexon=1.0, subexon_blocks=2, \
block_type='last_block', subexon_block_start=9, subexon_block_end=10, \
gap_block_start=3, gap_block_end=9)]
    """
    first_block = -1
    last_block = -1
    subexon_blocks = 0
    previous_subexon = np.nan
    problematic_list = []
    for (i, (subexon, _, _)) in enumerate(zip(subexons, starts, ends)):
        if not np.isnan(subexon):
            if not _same_subexon(subexon, previous_subexon):
                _store_problematic(problematic_list, sequence,
                                   previous_subexon, subexon_blocks, starts,
                                   ends, first_block, last_block,
                                   max_res_block, min_gap_block)
                previous_subexon = subexon
                first_block = i
                subexon_blocks = 1
            else:
                subexon_blocks += 1
                last_block = i
    _store_problematic(problematic_list, sequence, previous_subexon,
                       subexon_blocks, starts, ends, first_block, last_block,
                       max_res_block, min_gap_block)
    _choose_one(problematic_list)
    return problematic_list


def _complete_block_coords(block):
    if block.block_type == "first_block":
        return range(block.subexon_block_start, block.gap_block_end)
    return range(block.gap_block_start, block.subexon_block_end)


def cluster_subexon_blocks(blocks, margin=2):
    """
    Cluster sub-exon blocks starting or ending in the same column.
    """
    clusters = []
    columns = []
    for block in blocks:
        cluster_number = -1
        block_coords = _complete_block_coords(block)
        start = block_coords.start
        end = block_coords.stop
        for (index, cols) in enumerate(columns):
            if (start in cols) or (end in cols):
                cluster_number = index
        if cluster_number == -1:
            clusters.append([block])
            columns.append(range(max(0, start - margin), end + margin))
        else:
            clusters[cluster_number].append(block)
            cols = columns[cluster_number]
            columns[cluster_number] = range(
                max(0, min(cols.start, start - margin)),
                max(cols.stop, end + margin))
    return _merge_block_clusters(clusters, columns)


def _merge_block_clusters(clusters, columns):
    """
    Merge cluster sharing columns, return the list of clusters.
    """
    n_clusters = len(clusters)
    to_delete = set([])
    for actual_cluster in range(1, n_clusters):
        cols = set(columns[actual_cluster])
        for previous_cluster in range(actual_cluster):
            prev_cols = columns[previous_cluster]
            if previous_cluster not in to_delete:
                if cols.intersection(prev_cols):
                    clusters[previous_cluster].extend(clusters[actual_cluster])
                    columns[previous_cluster] = range(
                        max(0, min(min(cols), prev_cols.start)),
                        max(max(cols), prev_cols.end))
                    to_delete.update(actual_cluster)
                    break
    return [
        clus for (index, clus) in enumerate(clusters) if index not in to_delete
    ]


def problematic_block_clusters(
        msa_matrix,
        max_res_block=2,  # to move
        min_gap_block=2  # to consider
):
    """
    Return a Dict of problematic sub-exon block clusters.
    """
    blocks = []
    for (i, matrix_row) in enumerate(msa_matrix):
        subexons, starts, ends = resume_seq(matrix_row)
        blocks.extend(
            problematic_subexon_blocks(subexons,
                                       starts,
                                       ends,
                                       i,
                                       max_res_block=max_res_block,
                                       min_gap_block=min_gap_block))
    return cluster_subexon_blocks(blocks)


def move_block(  # pylint: disable=too-many-arguments
        matrix, row, block_start, block_end, destination_start,
        destination_end):
    """
    Move in-place a block in the (numpy) matrix.

    >>> import numpy as np
    >>> mat = np.matrix([[0,1,2,3,4], [5,6,7,8,9]])
    >>> move_block(mat, 0, 0, 2, 3, 5)
    matrix([[3, 4, 2, 0, 1],
            [5, 6, 7, 8, 9]])
    >>> move_block(mat, 1, 0, 2, 3, 5)
    matrix([[3, 4, 2, 0, 1],
            [8, 9, 7, 5, 6]])
    """
    subexon_block = np.copy(matrix[row, block_start:block_end])
    gap_block = np.copy(matrix[row, destination_start:destination_end])
    matrix[row, block_start:block_end] = gap_block
    matrix[row, destination_start:destination_end] = subexon_block
    return matrix


def move_subexon_block(msa_numpy, msa_matrix, subexon_block):
    """
    Swap a sub-exon block within its gap block.

    The movement depends on the problematic sub-exon block position, e.g.:

    ```
    from: [1.0, 1.0, 1.0, nan, nan, nan, 1.0, 2.0, 2.0, 2.0] # last_block
    to:   [1.0, 1.0, 1.0, 1.0, nan, nan, nan, 2.0, 2.0, 2.0]

    from: [1.0, 1.0, 1.0, 2.0, nan, nan, nan, 2.0, 2.0, 2.0] # first_block
    to:   [1.0, 1.0, 1.0, nan, nan, nan, 2.0, 2.0, 2.0, 2.0]
    ```
    """
    block_start = subexon_block.subexon_block_start
    block_end = subexon_block.subexon_block_end
    block_len = block_end - block_start
    if subexon_block.block_type == 'first_block':
        destination_start = subexon_block.gap_block_end - block_len
        destination_end = subexon_block.gap_block_end
    else:
        destination_start = subexon_block.gap_block_start
        destination_end = subexon_block.gap_block_start + block_len
    move_block(msa_numpy, subexon_block.sequence_index, block_start, block_end,
               destination_start, destination_end)
    move_block(msa_matrix, subexon_block.sequence_index, block_start,
               block_end, destination_start, destination_end)


def get_subexon_boundaries(col_clusters):
    """
    Returns a dict from (sequence number, subexon id) to start and end columns.
    """
    boundaries = dict()
    for cluster in col_clusters:
        for (seq, subexon) in enumerate(cluster.consensus):
            if not np.isnan(subexon):
                if (seq, subexon) in boundaries:
                    boundaries[(seq, subexon)][1] = cluster.end
                else:
                    boundaries[(seq, subexon)] = [cluster.start, cluster.end]
    return boundaries


def score_solution(msa_matrix):
    """
    Return the number of times a sub-exon cross an s-exon boundary.
    """
    col_patterns = column_patterns(msa_matrix)
    col_clusters = column_clusters(col_patterns)
    n_s_exons = len(col_clusters)
    n_sequences = len(col_clusters[0].consensus)
    score = 0
    if (n_s_exons > 1) and (n_sequences > 1):
        subexon_boundaries = get_subexon_boundaries(col_clusters)
        for cluster in range(1, n_s_exons):
            for coords in subexon_boundaries.values():
                start = col_clusters[cluster].start
                if coords[0] < start <= coords[1]:
                    score += 1
    return score


def move_problematic_block_clusters(
        msa,
        msa_matrix,
        max_res_block=2,  # to move
        min_gap_block=2  # to consider
):
    """
    Find and move all the problematic block clusters in a chimeric MSA.
    """
    block_clusters = problematic_block_clusters(msa_matrix,
                                                max_res_block=max_res_block,
                                                min_gap_block=min_gap_block)
    msa_numpy = np.array([list(rec) for rec in msa], np.character)
    for subexon_blocks in block_clusters:
        msa_copy = np.copy(msa_numpy)
        matrix_copy = np.copy(msa_matrix)
        starting_score = score_solution(matrix_copy)
        for subexon_block in subexon_blocks:
            move_subexon_block(msa_copy, matrix_copy, subexon_block)
        actual_score = score_solution(matrix_copy)
        if actual_score < starting_score:
            msa_numpy[:, :] = msa_copy
            msa_matrix[:, :] = matrix_copy
    return (msa_numpy, msa_matrix)


def _ask_directions(left_consensus, center_consensus, right_consensus):
    """
    Return a vector with the directions in which the residue should be moved.

    1: Move the residue to the gap column inserted at the right.
    0: Keep the residue/gap position.

    It returns an empty list if the residue can not be moved.

    >>> import numpy as np
    >>> _ask_directions([0.0, 2.0, 4.0], [0.0, 3.0, np.nan], [1.0, 3.0, 5.0])
    [0, 1, 0]
    >>> _ask_directions([0.0, 2.0, 4.0], [9.0, 3.0, np.nan], [1.0, 3.0, 5.0])
    []
    """
    directions = []
    for (left, center, right) in zip(left_consensus, center_consensus,
                                     right_consensus):
        direction = -1
        if np.isnan(center):
            direction = 0
        if center == right:
            direction = 1
        if center == left:
            direction = 0
        if direction == -1:
            return []
        directions.append(direction)
    return directions


def disintegration(msa_numpy, msa_matrix):
    """
    Disintegrate a one-residue length s-exons.

    It merges their residues to the contiguous s-exons.
    Return the modified MSA and sub-exon matrix.
    """
    col_patterns = column_patterns(msa_matrix)
    col_clusters = column_clusters(col_patterns)
    n_clusters = len(col_clusters)
    if n_clusters > 2:
        for index in range(1, n_clusters - 1):
            cluster = col_clusters[index]
            if cluster.start == cluster.end:
                left_cluster = col_clusters[index - 1]
                right_cluster = col_clusters[index + 1]
                directions = _ask_directions(left_cluster.consensus,
                                             cluster.consensus,
                                             right_cluster.consensus)
                if not directions:
                    continue
                col = cluster.start
                msa_matrix = np.insert(msa_matrix, col + 1, np.nan,
                                       axis=1)  # can not be in-place
                msa_numpy = np.insert(msa_numpy, col + 1, "-", axis=1)
                for (row, direction) in enumerate(directions):
                    if direction:  # is 1, then move
                        msa_matrix[row, col + 1] = msa_matrix[row, col]
                        msa_matrix[row, col] = np.nan
                        msa_numpy[row, col + 1] = msa_numpy[row, col]
                        msa_numpy[row, col] = "-"
                msa_numpy, msa_matrix = disintegration(msa_numpy, msa_matrix)
            else:
                continue
    return (msa_numpy, msa_matrix)
