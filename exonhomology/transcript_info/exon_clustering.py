from skbio.alignment import StripedSmithWaterman
from skbio.alignment._pairwise import blosum50


def exon_clustering(
        trx_data,  # pylint: disable=too-many-arguments
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

    nrows = trx_data.shape[0]

    trx_data = trx_data.assign(
        Cluster=0,
        QueryExon='',
        InputOrder=range(nrows),
        ProteinSequences=lambda df: df['Exon protein sequence'].map(
            lambda seq: str(seq).replace('*', '')
            ).values
    )
    # ?assign: For Python 3.5 and below, the order of keyword arguments is not
    # specified, you cannot refer to newly created or modified columns. All
    # items are computed first, and then assigned in alphabetical order.
    trx_data = trx_data.assign(
        SeqLength=lambda df: df['ProteinSequences'].map(len).values)

    # Sort by 'SeqLength'. We then sort by 'Exon stable ID' sequences with
    # equal length to ensure reproducibility:
    trx_data.sort_values(['SeqLength', 'Exon stable ID'],
                         inplace=True,
                         ascending=False)

    cluster_count = 0
    for i in range(nrows):
        i_index = trx_data.index[i]
        len_i = trx_data.at[i_index, 'SeqLength']
        if len_i < minimum_len:
            continue

        seq_i = trx_data.at[i_index, 'ProteinSequences']

        query_exon = trx_data.at[i_index, 'Exon stable ID']
        query = StripedSmithWaterman(
            seq_i,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
            substitution_matrix=substitution_matrix)

        cluster = trx_data.at[i_index, 'Cluster']
        if cluster == 0:
            cluster_count += 1
            trx_data.at[i_index, 'Cluster'] = cluster_count
            trx_data.at[i_index, 'QueryExon'] = query_exon
            cluster = cluster_count

        for j in range(i + 1, nrows):
            j_index = trx_data.index[j]
            len_j = trx_data.at[j_index, 'SeqLength']
            if len_j < minimum_len:
                continue
            seq_j = trx_data.at[j_index, 'ProteinSequences']

            aln = query(seq_j)

            aln_len = 0
            aln_len_i = 0
            aln_len_j = 0
            identical = 0
            for res_i, res_j in zip(aln.aligned_query_sequence,
                                    aln.aligned_target_sequence):
                aln_len += 1
                if res_i != '-':
                    aln_len_i += 1
                if res_j != '-':
                    aln_len_j += 1
                if res_i == res_j:
                    if res_i != '-':
                        identical += 1
                    else:
                        aln_len -= 1

            percent_identity = 100.0 * (identical / aln_len)
            # Because of the sort by SeqLength with ascending=False
            # j is always <= i, i.e. j is the shortest sequence of the pair
            coverage = 100.0 * (aln_len_j / len_j)

            if (coverage >= coverage_cutoff
                    and percent_identity >= percent_identity_cutoff):
                trx_data.at[j_index, 'Cluster'] = cluster
                trx_data.at[j_index, 'QueryExon'] = query_exon

    trx_data.sort_values('InputOrder', inplace=True)
    trx_data.drop(['InputOrder', 'ProteinSequences', 'SeqLength'],
                  axis=1,
                  inplace=True)

    return trx_data
