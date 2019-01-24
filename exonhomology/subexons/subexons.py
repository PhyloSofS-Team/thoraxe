"""
subexons: Module to identify and cut subexons.

This defines functions that take as input the DataFrame output from
transcript_info.read_transcript_info.
"""

import collections

import pandas as pd

from exonhomology import transcript_info
from exonhomology.transcript_info import phases

# A basic interval type, containing start and end coordinates
# and the set of components that are included on it.
Interval = collections.namedtuple('Interval', ['start', 'end', 'components'])


def _check_disjoint_intervals_inputs(starts, ends, ids):
    """Check the inputs of disjoint_intervals."""
    assert len(starts) == len(ends) == len(ids), 'starts, ends and ids '\
        'should have the same length.'
    assert all(starts[i] <= starts[i + 1] for i in range(len(starts) - 1)),\
        'starts should be sorted.'


def disjoint_intervals(starts, ends, ids):
    r"""
    Return list of disjoint intervals.

    Produce the set of disjoint intervals from a list of start and end
    positions (same chr) and optionally ids.
    Empty intervals are reported as well.
    !!! The intervals need to be sorted by start and end position before
    calling the function and the coordinates should be inclusive.

    >>> DisjointIntervals([1,10,20,40], [9,30,30,50], [0,1,2,3])
    [Interval(start=1, end=9, component=set([0])), \
    Interval(start=10, end=19, component=set([1])), \
    Interval(start=20, end=30, component=set([1, 2])), \
    Interval(start=31, end=39, component=set([])), \
    Interval(start=40, end=50, component=set([3]))]
    >>> DisjointIntervals([1,1], [9,7], [0,1])
    [Interval(start=1, end=7, component=set([0, 1])), \
    Interval(start=8, end=9, component=set([0]))]
    >>> DisjointIntervals([1,3], [7,9], [0,1])
    [Interval(start=1, end=2, component=set([0])), \
    Interval(start=3, end=7, component=set([0, 1])), \
    Interval(start=8, end=9, component=set([1]))]
    >>> DisjointIntervals([3,1], [9,7], [0,1])
    Traceback (most recent call last):
        ...
    AssertionError: starts should be sorted.
    >>> DisjointIntervals([3,1], [9,7,10], [0,1])
    Traceback (most recent call last):
        ...
    AssertionError: starts, ends and ids should have the same length.
    """
    _check_disjoint_intervals_inputs(starts, ends, ids)
    d_j = {}
    for i, start in enumerate(starts):
        d_j.setdefault(start, []).append((1, ids[i]))
    for i, end in enumerate(ends):
        # The end occurs after the end position,
        # explains also the BED style coding
        d_j.setdefault(end + 1, []).append((-1, ids[i]))
    j_pos = sorted(d_j.keys())
    interval_reduce = []
    # The list of open intervals from the last turn
    last = set([])
    for left, right in zip(j_pos[:-1], j_pos[1:]):
        cint = Interval(left, right - 1, set(last))
        cint.components.update([x[1] for x in d_j[left] if x[0] == 1])
        cint.components.symmetric_difference_update(
            [x[1] for x in d_j[left] if x[0] == -1])
        interval_reduce.append(cint)
        last = cint.components
    # Need to decide of the format.
    return interval_reduce


def _get_subexon_seq(exon_row, interval):
    """Return the interval sequence of the subexon."""
    seq = exon_row['Exon sequence']
    if exon_row['Strand'] == 1:
        exon_start = exon_row['Exon region start (bp)']
        interval_start = interval.start - exon_start
        # Because interval.end includes the last base (close interval)
        # and Python indexing doesn't include the last value, +1 is needed:
        interval_end = interval.end - exon_start + 1
    else:
        exon_start = exon_row['Exon region end (bp)']
        interval_start = exon_start - interval.end
        interval_end = exon_start - interval.start + 1

    assert interval_start < interval_end, 'Interval end should not be '\
        'lower than its start'

    return seq[interval_start:interval_end]


def _subexon_ranks(strand, transcript_len):
    """
    Return a list of the subexon ranks.

    NOTE: Rank starts in 0 to make this compatible with
    end_phase_previous_exon that expect exon_pos(ition) to start in 0.

    >>> _subexon_ranks(-1, 10)
    [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
    >>> _subexon_ranks(1, 10)
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    """
    if strand == -1:
        return list(range(transcript_len - 1, -1, -1))

    return list(range(0, transcript_len))


def _sort_subexons_coordinates(subexon_df):
    """Sort subexons and return the DataFrame."""
    transcript_data_frames = []

    id_columns = ['Gene stable ID', 'Transcript stable ID']
    for _, transcript_df in subexon_df.groupby(id_columns):
        transcript_data_frames.append(
            transcript_df.sort_values(
                by='Genomic coding start', ascending=True))

    return pd.concat(transcript_data_frames)


def _add_subexon_rank(subexon_df):
    """
    Sort subexon_df and add a column with the subexon rank in the transcript.

    It returns the modified pandas' DataFrame.
    """
    subexon_df = _sort_subexons_coordinates(subexon_df)

    ranks = []
    id_columns = ['Gene stable ID', 'Transcript stable ID']
    for _, transcript_df in subexon_df.groupby(id_columns):
        strand = transcript_df.loc[transcript_df.index[0], 'Strand']
        transcript_len = transcript_df.shape[0]
        ranks.extend(_subexon_ranks(strand, transcript_len))

    subexon_df['Subexon rank in transcript'] = ranks

    return subexon_df


def _subexon_info(gene_df, exon_df, intervals):
    """Return a pandas' DataFrame for subexons."""
    # We are going to use transcript_info functions to transcribe and
    # merge subexons. Therefore, we are going to use the exon table column
    # names (column_names keys) at the beginning, and then we are going to
    # change this column names by the column_names values before returning
    # subexon_df, columns with '' as name won't be in the output :
    column_names = {
        'Subexon ID': 'Subexon ID',
        'Interval number': 'Interval number',
        'Exon sequence': 'Subexon sequence',
        'Genomic coding start': 'Subexon genomic coding start',
        'Genomic coding end': 'Subexon genomic coding end',
        # Transcript information
        'Gene stable ID': 'Gene stable ID',
        'Transcript stable ID': 'Transcript stable ID',
        'Strand': 'Strand',
        'Exon rank in transcript': 'Exon rank in transcript',
        'Transcript stable ID cluster': 'Transcript stable ID cluster',
        # 'Exon stable ID' is kept to merge the exon and subexon tables later
        'Exon stable ID': 'Exon stable ID',
        'Exon stable ID cluster': 'Exon stable ID cluster',
        'Cluster': 'Cluster'
    }
    # cDNA coding start and end are only needed for _get_exon_cds
    column_names_to_delete = ['cDNA coding start', 'cDNA coding end']

    subexon_df = pd.DataFrame(
        columns=list(column_names.keys()) + column_names_to_delete)

    for interval_number, interval in enumerate(intervals):
        for exon_id in interval.components:
            exon_row = exon_df.loc[exon_id]
            seq = _get_subexon_seq(exon_row, interval)
            subexon_info = {
                'Subexon ID': exon_id + '_SE_' + str(interval_number),
                'Interval number': interval_number,
                'Exon stable ID': exon_id,
                'Genomic coding start': interval.start,
                'Genomic coding end': interval.end,
                'Cluster': exon_row['Cluster'],

                # CDS sequence of the subexon:
                'Exon sequence': seq,

                # cDNA coding start and end have dummy values to make
                # _get_exon_cds to return the input sequence without
                # modifications:
                'cDNA coding start': 0,
                'cDNA coding end': len(seq) - 1
            }
            for _, exon_row in gene_df[gene_df['Exon stable ID'] ==
                                       exon_id].iterrows():
                transcript_row_info = {
                    'Gene stable ID':
                    exon_row['Gene stable ID'],
                    'Transcript stable ID':
                    exon_row['Transcript stable ID'],
                    'Transcript stable ID cluster':
                    exon_row['Transcript stable ID cluster'],
                    'Exon stable ID cluster':
                    exon_row['Exon stable ID cluster'],
                    'Strand':
                    exon_row['Strand'],
                    'Exon rank in transcript':
                    exon_row['Exon rank in transcript']
                }
                transcript_row_info.update(subexon_info)
                subexon_df = subexon_df.append(
                    transcript_row_info, ignore_index=True)

    subexon_df = _add_subexon_rank(subexon_df)
    subexon_df = subexon_df.sort_values(by=[
        'Gene stable ID', 'Transcript stable ID', 'Subexon rank in transcript'
    ])

    _add_phases(
        subexon_df,
        rank_column='Subexon rank in transcript',
        seq_column='Exon sequence')

    transcript_info.add_protein_seq(subexon_df)

    transcript_info.merge_identical_exons(
        subexon_df, exon_id_column='Subexon ID')

    subexon_df.drop(columns=column_names_to_delete, inplace=True)
    subexon_df.rename(columns=column_names, inplace=True)

    return subexon_df


def _add_phases(subexon_df,
                rank_column='Subexon rank in transcript',
                seq_column='Subexon sequence'):
    """
    Add columns with start and end phases.

    This function need to be run after _add_subexon_rank is used to
    assign rank to the subexons in the transcript.

    subexon_df needs to be sorted by subexon rank.
    """
    subexon_df['Start phase'] = -99
    subexon_df['End phase'] = -99

    n_rows = subexon_df.shape[0]
    row_number = 0
    while row_number < n_rows:
        row_index = subexon_df.index[row_number]
        row = subexon_df.loc[row_index, :]

        cdna_len = len(row[seq_column])

        if row_number == 0:
            previous_end_phase = 0
        else:
            subexon_rank = row[rank_column]
            previous_end_phase = phases.end_phase_previous_exon(
                subexon_df, subexon_rank, subexon_df.index[row_number - 1])

        start_phase, end_phase = phases.calculate_phase(
            cdna_len, previous_end_phase)

        subexon_df.at[row_index, 'Start phase'] = start_phase
        subexon_df.at[row_index, 'End phase'] = end_phase

        row_number += 1


def _add_transcript_fraction(subexon_df):
    """
    Add the fraction of the transcripts in which the subexon is present.

    1.0 for constitutive exons, < 1.0 for alternative exons
    """
    subexon_counts = subexon_df.groupby([
        'Gene stable ID', 'Subexon ID'
    ]).size().to_frame('Number of transcripts for subexon').reset_index()

    transcript_counts = subexon_df.loc[:, [
        'Gene stable ID', 'Transcript stable ID'
    ]].drop_duplicates().groupby([
        'Gene stable ID'
    ]).size().to_frame('Transcripts in gene').reset_index()

    counts = pd.merge(
        subexon_counts, transcript_counts, how='inner', on='Gene stable ID')

    counts['Transcript fraction'] = (
        counts['Number of transcripts for subexon'] /
        counts['Transcripts in gene'].astype('float64'))

    return pd.merge(
        subexon_df, counts, how='inner', on=['Gene stable ID', 'Subexon ID'])


def _merge_intervals_with_one_element(intervals):
    """
    Merge contigous intervals with the same and unique components.

    Intervals should be sorted.
    """
    merged_intervals = []
    n_intervals = len(intervals)
    if n_intervals > 1:
        for index in range(n_intervals):
            if index == 0:
                merged_intervals.append(intervals[index])
            else:
                previous = merged_intervals[-1]
                actual = intervals[index]
                print(previous)
                print(actual)
                if len(actual.components) == 1 and len(
                        previous.components) == 1 and sorted(
                            actual.components) == sorted(
                                previous.
                                components) and actual.start - actual.end == 1:
                    merged_intervals[-1] = Interval(previous.start, actual.end,
                                                    actual.components)
                else:
                    merged_intervals.append(actual)
    else:
        return intervals
    return merged_intervals


def create_subexon_table(transcript_info):
    """Return a subexon table."""
    subexon_data_frames = []
    for _, gene_df in transcript_info.groupby('Gene stable ID'):

        # We need the genomic coding coordinates and the sequence of each exon
        exon_df = gene_df.drop_duplicates(subset='Exon stable ID')
        exon_df.set_index('Exon stable ID', inplace=True)
        exon_df = exon_df.sort_values(
            by=['Genomic coding start', 'Genomic coding end'])

        intervals = disjoint_intervals(
            list(exon_df['Genomic coding start']),
            list(exon_df['Genomic coding end']), list(exon_df.index))

        intervals = [
            interval for interval in intervals if len(interval.components) > 0
        ]

        intervals = sorted(intervals, key=lambda interval: interval.start)

        if gene_df['Gene stable ID'].unique()[0] == 'ENSMMUG00000004060':
            print(intervals)

        intervals = _merge_intervals_with_one_element(intervals)

        if gene_df['Gene stable ID'].unique()[0] == 'ENSMMUG00000004060':
            print(intervals)

        subexon_df = _subexon_info(gene_df, exon_df, intervals)
        subexon_data_frames.append(subexon_df)

    subexon_table = pd.concat(subexon_data_frames)
    subexon_table = _add_transcript_fraction(subexon_table)

    return subexon_table
