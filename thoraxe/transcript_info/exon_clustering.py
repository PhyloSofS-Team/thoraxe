"""Functions to cluster exons."""

import warnings

import numpy as np

from skbio.alignment import StripedSmithWaterman
from skbio.alignment._pairwise import blosum50


def coverage(seq, seq_len):
    """
    Return coverage of the sequence in the alignment.

    >>> coverage("AAAA----", 8)
    50.0
    """
    return 100.0 * (len(seq) - seq.count('-')) / seq_len


def percent_identity(query, target):
    """
    Return percent identity of the aligned sequnces.

    >>> percent_identity("AA---", "AAAA-")
    50.0
    """
    aln_len = len(query)
    if aln_len != len(target):
        warnings.warn('{} and {} have different lengths'.format(query, target))
        return np.nan

    identical = 0
    for res_i, res_j in zip(query, target):
        if res_i == res_j:
            if res_i != '-':
                identical += 1
            else:
                aln_len -= 1

    return 100.0 * (identical / aln_len)


def exon_clustering(  # pylint: disable=too-many-arguments,too-many-locals
        trx_data,
        minimum_len=4,
        coverage_cutoff=80.0,
        percent_identity_cutoff=30.0,
        gap_open_penalty=10,
        gap_extend_penalty=1,
        substitution_matrix=None):
    """
    Cluster exons based on their sequence identity after local alignment.

    It uses a Hobohm I sequence clustering algorithm to perform a fast
    clustering of the exons. Exons are sorted from the longest to the
    shortest before start the clustering, but at the end the returned
    table as the same order as the input. The returned data frame as two
    extra columns, Cluster and QueryExon that contains the cluster number
    and the name of the query sequence used to join that one into that
    cluster.

    The alignment is performed using StripedSmithWaterman from scikit-bio.
    The keyword arguments gap_open_penalty, gap_extend_penalty and
    substitution_matrix are passed to StripedSmithWaterman.

    Exons with a length less than the minimum_len (default: 4) are not
    clustered. Non-clustered exons have Cluster number 0 and an empty
    string as ExonQuery.

    The coverage and percent identity cutoff to decide if an exon sequence
    belongs to one cluster can be modified with the keyword arguments
    coverage_cutoff and percent_identity_cutoff (default to 80.0 and
    30.0, respectively).
    """
    if substitution_matrix is None:
        substitution_matrix = blosum50

    nrows = len(trx_data)

    trx_data = trx_data.assign(
        Cluster=0,
        QueryExon='',
        TargetCoverage=np.nan,
        PercentIdentity=np.nan,
        AlignedQuery='',
        AlignedTarget='',
        InputOrder=range(nrows),
        ProteinSequences=lambda df: df['ExonProteinSequence'].map(
            lambda seq: str(seq).replace('*', '')).values)
    # ?assign: For Python 3.5 and below, the order of keyword arguments is not
    # specified, you cannot refer to newly created or modified columns. All
    # items are computed first, and then assigned in alphabetical order.
    trx_data = trx_data.assign(
        SeqLength=lambda df: df['ProteinSequences'].map(len).values)

    # Sort by 'SeqLength'. We then sort by 'ExonStableID' sequences with
    # equal length to ensure reproducibility:
    trx_data.sort_values(['SeqLength', 'ExonStableID'],
                         inplace=True,
                         ascending=False)

    row_list = trx_data.to_dict('records')

    cluster_count = 0
    for i in range(nrows):
        if row_list[i]['SeqLength'] < minimum_len:
            continue
        query = StripedSmithWaterman(
            row_list[i]['ProteinSequences'],
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
            substitution_matrix=substitution_matrix)

        query_exon = row_list[i]['ExonStableID']
        i_index = trx_data.index[i]
        cluster = trx_data.at[i_index, 'Cluster']
        if cluster == 0:
            cluster_count += 1
            trx_data.at[i_index, 'Cluster'] = cluster_count
            trx_data.at[i_index, 'QueryExon'] = query_exon
            cluster = cluster_count

        for j in range(i + 1, nrows):
            j_index = trx_data.index[j]
            if row_list[j]['SeqLength'] < minimum_len or trx_data.at[
                    j_index, 'Cluster'] != 0:
                continue

            aln = query(row_list[j]['ProteinSequences'])

            # Because of the sort by SeqLength with ascending=False, the target
            # j is always <= i, i.e. j is the shortest sequence of the pair.
            target_coverage = coverage(aln.aligned_target_sequence,
                                       row_list[j]['SeqLength'])

            pid = percent_identity(aln.aligned_query_sequence,
                                   aln.aligned_target_sequence)

            if (target_coverage >= coverage_cutoff
                    and pid >= percent_identity_cutoff):
                trx_data.at[j_index, 'Cluster'] = cluster
                trx_data.at[j_index, 'QueryExon'] = query_exon
                trx_data.at[j_index, 'TargetCoverage'] = target_coverage
                trx_data.at[j_index, 'PercentIdentity'] = pid
                trx_data.at[j_index,
                            'AlignedQuery'] = aln.aligned_query_sequence
                trx_data.at[j_index,
                            'AlignedTarget'] = aln.aligned_target_sequence

    trx_data.sort_values('InputOrder', inplace=True)
    trx_data.drop(['InputOrder', 'ProteinSequences', 'SeqLength'],
                  axis=1,
                  inplace=True)

    return trx_data
