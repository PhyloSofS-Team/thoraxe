"""
Functions for the rescue phase.
"""

import collections

import pandas as pd
import numpy as np
from skbio.alignment import StripedSmithWaterman
from skbio.alignment._pairwise import blosum50
from thoraxe import transcript_info


def _get_subexons_to_rescue(subexon_table):
    """Return a DataFrame with the subexons to rescue."""
    return subexon_table.loc[subexon_table['Cluster'] < 0, [
        'Subexon ID cluster', 'Cluster', 'Exon protein sequence'
    ]].drop_duplicates(subset='Subexon ID cluster').to_dict('records')


def _get_sequence_list(sequences, minimum_len):
    """Helper funtion to return a set of unique sequences."""
    seq_list = set([])
    for seq in sequences:
        seq = str(seq).replace('*', '')
        if len(seq) > minimum_len:
            seq_list.add(seq)
    return seq_list


def _get_cluster2sequence(subexon_table, minimum_len=0):
    """Return a dict from cluster id to cluster unique sequences."""
    clusterized = subexon_table.loc[subexon_table['Cluster'] > 0, [
        'Subexon ID cluster', 'Cluster', 'Exon protein sequence'
    ]]
    clusterized.drop_duplicates('Subexon ID cluster', inplace=True)
    clusterized.drop(columns='Subexon ID cluster', inplace=True)
    cluster2sequence = clusterized.groupby('Cluster').agg(
        lambda col: _get_sequence_list(col, minimum_len)).to_dict('index')
    return {
        cluster: seqs['Exon protein sequence']
        for cluster, seqs in cluster2sequence.items()
        if seqs['Exon protein sequence']
    }


def _get_cluster_number(aln_stats):
    """Helper function to get the cluster id with the best matches."""
    df_stats = pd.DataFrame(aln_stats)
    df_stats['number'] = 1
    df_stats = df_stats.groupby('cluster').agg({
        'coverage': max,
        'percent identity': max,
        'number': sum
    })
    df_stats.sort_values(
        by=['number', 'coverage', 'percent identity'],
        ascending=False,
        inplace=True)
    return df_stats.index[0]


def _get_aln_stats(  # pylint: disable=too-many-arguments
        query, query_len, origin_cluster, cluster2sequence, coverage_cutoff,
        percent_identity_cutoff):
    """Helper function to get a dict with the aln stats of the matches."""
    aln_stats = collections.defaultdict(list)
    for cluster, sequences in cluster2sequence.items():
        if cluster == origin_cluster:
            continue
        for sequence in sequences:
            aln = query(sequence)
            pid = transcript_info.percent_identity(aln.aligned_query_sequence,
                                                   aln.aligned_target_sequence)
            cov = transcript_info.coverage(aln.aligned_query_sequence,
                                           query_len)
            if pid >= percent_identity_cutoff and cov >= coverage_cutoff:
                aln_stats['percent identity'].append(pid)
                aln_stats['coverage'].append(cov)
                aln_stats['cluster'].append(cluster)
    return aln_stats


def _get_subexon2cluster(  # pylint: disable=too-many-arguments
        to_rescue,
        cluster2sequence,
        coverage_cutoff=80.0,
        percent_identity_cutoff=30.0,
        gap_open_penalty=10,
        gap_extend_penalty=1,
        substitution_matrix=None):
    """
    Return a dict from 'Subexon ID cluster' to cluster.

    This function takes
    """
    if substitution_matrix is None:
        substitution_matrix = blosum50

    subexon2cluster = {}
    for row in to_rescue:
        origin_cluster = -1 * row['Cluster']
        seq = str(row['Exon protein sequence']).replace('*', '')
        query = StripedSmithWaterman(
            seq,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
            substitution_matrix=substitution_matrix)
        aln_stats = _get_aln_stats(query, len(seq), origin_cluster,
                                   cluster2sequence, coverage_cutoff,
                                   percent_identity_cutoff)
        if aln_stats:
            subexon2cluster[row['Subexon ID cluster']] = _get_cluster_number(
                aln_stats)

    return subexon2cluster


def modify_subexon_cluster(subexon_table, subexon2cluster):
    """Modify subexon_table using the mapping in subexon2cluster."""
    for subexon, cluster in subexon2cluster.items():
        subexon_index = subexon_table.index[subexon_table['Subexon ID cluster']
                                            == subexon]
        for i in subexon_index:
            subexon_table.at[i, 'Cluster'] = cluster
            subexon_table.at[i, 'QueryExon'] = ''
            subexon_table.at[i, 'AlignedQuery'] = ''
            subexon_table.at[i, 'AlignedTarget'] = ''
            subexon_table.at[i, 'TargetCoverage'] = np.nan
            subexon_table.at[i, 'PercentIdentity'] = np.nan

    return subexon_table


def subexon_rescue_phase(  # pylint: disable=too-many-arguments
        subexon_table,
        minimum_len=4,
        coverage_cutoff=80.0,
        percent_identity_cutoff=30.0,
        gap_open_penalty=10,
        gap_extend_penalty=1,
        substitution_matrix=None):
    """
    Execute the subexon rescue phase.

    This function looks the subexons to rescues in the subexon_table by
    looking at the 'Cluster' numbers that are lower than 0, i.e. -1 means that
    the subexon has been cleaned up from 'Cluster' 1.
    Then, each subexon to rescue is aligned against the subexons in other
    'Cluster's, different from the original one. If it is nicely aligned
    against another cluster, the 'Cluster' number is modified in the
    subexon_table to match the new proposed 'Cluster' and the data relative to
    the original pairwise alignment are cleaned up.
    This function returns the list of 'Cluster's that have been modified.
    """
    to_rescue = _get_subexons_to_rescue(subexon_table)
    cluster2sequence = _get_cluster2sequence(subexon_table, minimum_len)
    subexon2cluster = _get_subexon2cluster(
        to_rescue,
        cluster2sequence,
        coverage_cutoff=coverage_cutoff,
        percent_identity_cutoff=percent_identity_cutoff,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty,
        substitution_matrix=substitution_matrix)
    modify_subexon_cluster(subexon_table, subexon2cluster)
    return list(subexon2cluster.values())
