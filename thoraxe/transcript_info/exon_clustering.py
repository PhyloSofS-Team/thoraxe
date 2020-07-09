"""
Exon clustering: Functions to cluster exons using pairwise alignments.
"""

import collections
import warnings

import numpy as np

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum50

from thoraxe import subexons


def coverage(seq, seq_len):
    """
    Return coverage of the sequence in the alignment.

    >>> coverage("AAAA----", 8)
    50.0
    """
    res = len(seq.replace('-', ''))
    return 100.0 * res / seq_len


def coverage_shortest(seq_query, seq_target, seq_len):
    """
    Return coverage of the shortest sequence in the alignment.

    >>> coverage_shortest("AAAA----", "AAAAAAAA", 8)
    50.0
    """
    res_query = len(seq_query.replace('-', ''))
    res_target = len(seq_target.replace('-', ''))
    return 100.0 * min(res_query, res_target) / seq_len


def percent_identity(query, target):
    """
    Return percent identity of the aligned sequences.

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


def _align(seq_a, seq_b, gap_open_penalty, gap_extend_penalty,
           substitution_matrix):
    "Align the sequnces using the BioPython pairwise2 aligner."
    alignments = pairwise2.align.localds(seq_a,
                                         seq_b,
                                         substitution_matrix,
                                         gap_open_penalty,
                                         gap_extend_penalty,
                                         one_alignment_only=True)
    if alignments:
        aln = alignments[0]
        return aln[0], aln[1]
    return "", ""


def _max_percent_identity(residues_count_query, residues_count_target,
                          len_target):
    """
    Return the maximum percent identity between the two sequences.
    """
    if residues_count_query.keys().isdisjoint(residues_count_target.keys()):
        return 0
    residues_count = residues_count_target.copy()
    residues_count.subtract(residues_count_query)
    return 100.0 * (len_target - sum(
        val for val in residues_count.values() if val > 0)) / len_target


def _merge_clusters(trx_data, clusters_to_merge):
    """
    Merge the clusters in the transcript table.
    """
    old2new = dict()
    for group in clusters_to_merge:
        clusters = sorted(group)
        n_clusters = len(clusters)
        if n_clusters > 1:
            new_cluster_id = clusters[0]
            for i in range(1, n_clusters):
                old2new[clusters[i]] = new_cluster_id
    nrows = len(trx_data)
    for i in range(nrows):
        i_index = trx_data.index[i]
        cluster = trx_data.at[i_index, 'Cluster']
        trx_data.at[i_index, 'Cluster'] = old2new.get(cluster, cluster)
    return trx_data


def exon_clustering(  # noqa pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
        trx_data,
        minimum_len=4,
        coverage_cutoff=80.0,
        percent_identity_cutoff=30.0,
        gap_open_penalty=-10,
        gap_extend_penalty=-1,
        substitution_matrix=None,
        merge_clusters=False):
    """
    Cluster exons based on their sequence identity after local alignment.

    It uses a Hobohm I sequence clustering algorithm to perform a fast
    clustering of the exons. Exons are sorted from the longest to the
    shortest before start the clustering, but at the end the returned
    table as the same order as the input. The returned data frame as two
    extra columns, Cluster and QueryExon that contains the cluster number
    and the name of the query sequence used to join that one into that
    cluster.

    The alignment is performed using Smith-Waterman from BioPython.
    The keyword arguments gap_open_penalty, gap_extend_penalty and
    substitution_matrix are passed to pairwise2.

    Exons with a length less than the minimum_len (default: 4) are not
    clustered. Non-clustered exons have Cluster number 0 and an empty
    string as ExonQuery.

    The coverage and percent identity cutoff to decide if an exon sequence
    belongs to one cluster can be modified with the keyword arguments
    coverage_cutoff and percent_identity_cutoff (default to 80.0 and
    30.0, respectively).

    Exons are aligned against others until they reach a percent identity
    higher than percent_identity_cutoff plus 30, to allow low identity matches
    to be replaced by better ones.
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

    # Sort by 'SeqLength'. We then sort by 'ExonID' sequences with
    # equal length to ensure reproducibility:
    trx_data.sort_values(['SeqLength', 'ExonID'],
                         inplace=True,
                         ascending=False)

    trx_data['KeepSearching'] = 0

    row_list = trx_data.to_dict('records')

    clusters_to_merge = []

    cluster_count = 0
    aln_memo = dict()
    residues_counts = [
        collections.Counter(row_list[i]['ProteinSequences'])
        for i in range(nrows)
    ]
    for i in range(nrows):
        len_i = row_list[i]['SeqLength']
        if len_i < minimum_len:
            continue

        query_exon = row_list[i]['ExonID']
        i_index = trx_data.index[i]
        cluster = trx_data.at[i_index, 'Cluster']
        if cluster == 0:
            cluster_count += 1
            trx_data.at[i_index, 'Cluster'] = cluster_count
            trx_data.at[i_index, 'QueryExon'] = query_exon
            cluster = cluster_count

        for j in range(i + 1, nrows):
            j_index = trx_data.index[j]
            len_j = row_list[j]['SeqLength']
            if ((len_j < minimum_len)
                    or ((trx_data.at[j_index, 'Cluster'] != 0)
                        and not trx_data.at[j_index, 'KeepSearching'])):
                continue

            target_exon = row_list[j]['ExonID']
            exon_pair_id = tuple(sorted([query_exon, target_exon]))
            if exon_pair_id in aln_memo:
                pid, target_coverage = aln_memo[exon_pair_id]
            else:
                max_pid = _max_percent_identity(residues_counts[i],
                                                residues_counts[j], len_j)
                if max_pid < percent_identity_cutoff:
                    aln_memo[exon_pair_id] = (0.0, 0.0)
                    continue
                aln_query, aln_target = _align(row_list[i]['ProteinSequences'],
                                               row_list[j]['ProteinSequences'],
                                               gap_open_penalty,
                                               gap_extend_penalty,
                                               substitution_matrix)
                target_coverage = coverage_shortest(aln_query, aln_target,
                                                    len_j)
                pid = percent_identity(aln_query, aln_target)
                aln_memo[exon_pair_id] = (pid, target_coverage)

            if (target_coverage >= coverage_cutoff
                    and pid >= percent_identity_cutoff):
                if (trx_data.at[j_index, 'KeepSearching']
                        and trx_data.at[j_index, 'PercentIdentity'] > pid):
                    continue
                if merge_clusters:
                    current_cluster = trx_data.at[j_index, 'Cluster']
                    if current_cluster != 0:
                        subexons.update_to_merge_list(clusters_to_merge,
                                                      current_cluster, cluster)
                trx_data.at[j_index, 'Cluster'] = cluster
                trx_data.at[j_index, 'QueryExon'] = query_exon
                trx_data.at[j_index, 'TargetCoverage'] = target_coverage
                trx_data.at[j_index, 'PercentIdentity'] = pid
                trx_data.at[j_index, 'AlignedQuery'] = aln_query
                trx_data.at[j_index, 'AlignedTarget'] = aln_target
                if pid <= min(100, percent_identity_cutoff + 30):
                    trx_data.at[j_index, 'KeepSearching'] = 1
                else:
                    trx_data.at[j_index, 'KeepSearching'] = 0

    if clusters_to_merge:
        trx_data = _merge_clusters(trx_data, clusters_to_merge)

    trx_data.sort_values('InputOrder', inplace=True)
    trx_data.drop(
        ['InputOrder', 'ProteinSequences', 'SeqLength', 'KeepSearching'],
        axis=1,
        inplace=True)

    return trx_data
