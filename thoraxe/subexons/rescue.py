"""
Functions for the rescue phase.
"""

import collections
import logging

import pandas as pd
import numpy as np
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum50

from thoraxe import transcript_info


def _get_subexons_to_rescue(subexon_table):
    """Return a DataFrame with the subexons to rescue."""
    return subexon_table.loc[
        subexon_table['Cluster'] < 0,
        ['GeneID', 'SubexonIDCluster', 'Cluster', 'SubexonProteinSequence'
         ]].drop_duplicates(subset='SubexonIDCluster').to_dict('records')


def _get_cluster_number(aln_stats):
    """Helper function to get the cluster id with the best matches."""
    df_stats = pd.DataFrame(aln_stats)
    df_stats['number'] = 1
    df_stats = df_stats.groupby('cluster').agg({
        'coverage': max,
        'percent identity': max,
        'number': sum
    })
    df_stats.sort_values(by=['number', 'coverage', 'percent identity'],
                         ascending=False,
                         inplace=True)
    return df_stats.index[0]


def get_unaligned_regions(seq_i, seq_j, minimum_len=4):
    """
    Return the ungapped regions in `seq_j` that do not align with `seq_i`.

    Returned regions should have at least `minimum_len`.

    >>> seq_i = '------MHKCLVDE------YTEDQGGFRK------'
    >>> seq_j = 'MLLHYHHHKC-------LMCYTRDLHG---IH-L-K'
    >>> get_unaligned_regions(seq_i, seq_j)
    ['MLLHYH', 'IHLK']
    >>> seq_i = 'XXXXXXXXXXMHKCLVDE------YTEDQGGFRK'
    >>> seq_j = '----MLLHYHHHKC-------LMCYTRDLHG---'
    >>> get_unaligned_regions(seq_i, seq_j)
    ['MLLHYH']
    >>> seq_i = 'XXXXXXXXXXMHKCLVDE------YTEDQGGFRK'
    >>> seq_j = '----MLLHYHHHKC-----XXLMCYTRDLHG---'
    >>> get_unaligned_regions(seq_i, seq_j)
    ['MLLHYH', 'XXLMC']
    """
    j_sequences = []
    j_region = []
    in_region = False
    for (res_i, res_j) in zip(seq_i, seq_j):
        if res_i in {'-', 'X'}:
            in_region = True
            if res_j != '-':
                j_region.append(res_j)
        elif in_region:
            seq = ''.join(j_region)
            j_region.clear()
            in_region = False
            if len(seq) >= minimum_len:
                j_sequences.append(seq)

    seq = ''.join(j_region)
    if len(seq) >= minimum_len:
        j_sequences.append(seq)

    return j_sequences


def _get_gene_cluster_to_unaligned(gene_id, msa, minimum_len):
    """
    Dict from (gene, cluster) to the sequences that doesn't align to the gene.
    """
    unaligned_seqs = set()
    n_seq = len(msa)
    for i in range(0, n_seq):
        seq_i = msa[i, :]
        if seq_i.id == gene_id:
            for j in range(0, n_seq):
                if i != j:
                    seq_j = msa[j, :]
                    unaligned_seqs.update(
                        get_unaligned_regions(seq_i,
                                              seq_j,
                                              minimum_len=minimum_len))
    return list(unaligned_seqs)


def _get_aln_stats(  # pylint: disable=too-many-arguments,too-many-locals
        row, cluster2data, coverage_cutoff, percent_identity_cutoff,
        gap_open_penalty, gap_extend_penalty, substitution_matrix,
        minimum_len):
    """Helper function to get a dict with the aln stats of the matches."""

    query_seq = str(row['SubexonProteinSequence']).replace('*', '')
    query_len = len(query_seq)

    aln_stats = collections.defaultdict(list)

    cluster2sequences = collections.defaultdict(list)

    origin_cluster = -1 * row['Cluster']
    query_gene = row['GeneID']
    for (cluster, (_, chimerics, msa)) in cluster2data.items():
        if cluster != origin_cluster:
            if query_gene not in chimerics:
                cluster2sequences[cluster] = [
                    seq.replace('-', '') for (seq, _) in chimerics.values()
                ]
            elif len(chimerics) > 1:
                cluster2sequences[cluster] = _get_gene_cluster_to_unaligned(
                    query_gene, msa, minimum_len)

    for cluster, sequences in cluster2sequences.items():
        if cluster == origin_cluster:
            continue
        for sequence in sequences:
            alignments = pairwise2.align.localds(query_seq,
                                                 sequence,
                                                 substitution_matrix,
                                                 gap_open_penalty,
                                                 gap_extend_penalty,
                                                 one_alignment_only=True)
            if alignments:
                aln = alignments[0]
                pid = transcript_info.percent_identity(aln[0], aln[1])
                cov = transcript_info.coverage(aln[0], query_len)
                if pid >= percent_identity_cutoff and cov >= coverage_cutoff:
                    aln_stats['percent identity'].append(pid)
                    aln_stats['coverage'].append(cov)
                    aln_stats['cluster'].append(cluster)
    return aln_stats


def _get_subexon2cluster(  # pylint: disable=too-many-arguments
        to_rescue,
        cluster2data,
        coverage_cutoff=80.0,
        percent_identity_cutoff=30.0,
        gap_open_penalty=-10,
        gap_extend_penalty=-1,
        substitution_matrix=None,
        minimum_len=4):
    """
    Return a dict from 'SubexonIDCluster' to cluster.

    This function takes
    """
    if substitution_matrix is None:
        substitution_matrix = blosum50

    subexon2cluster = {}
    for row in to_rescue:
        aln_stats = _get_aln_stats(row,
                                   cluster2data,
                                   coverage_cutoff,
                                   percent_identity_cutoff,
                                   gap_open_penalty=gap_open_penalty,
                                   gap_extend_penalty=gap_extend_penalty,
                                   substitution_matrix=substitution_matrix,
                                   minimum_len=minimum_len)
        if aln_stats:
            logging.warning(
                'Sub-exon %s from gene %s is being rescued (previous cluster: %s)',
                row["SubexonIDCluster"], row["GeneID"], -1 * row["Cluster"])
            subexon2cluster[row['SubexonIDCluster']] = _get_cluster_number(
                aln_stats)

    return subexon2cluster


def modify_subexon_cluster(subexon_table, subexon2cluster):
    """Modify subexon_table using the mapping in subexon2cluster."""
    for subexon, cluster in subexon2cluster.items():
        subexon_index = subexon_table.index[subexon_table['SubexonIDCluster']
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
        cluster2data,
        subexon_table,
        minimum_len=4,
        coverage_cutoff=80.0,
        percent_identity_cutoff=30.0,
        gap_open_penalty=-10,
        gap_extend_penalty=-1,
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
    subexon2cluster = _get_subexon2cluster(
        to_rescue,
        cluster2data,
        coverage_cutoff=coverage_cutoff,
        percent_identity_cutoff=percent_identity_cutoff,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty,
        substitution_matrix=substitution_matrix,
        minimum_len=minimum_len)
    modify_subexon_cluster(subexon_table, subexon2cluster)
    return list(subexon2cluster.values())
