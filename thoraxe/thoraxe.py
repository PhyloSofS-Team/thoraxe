"""
ThorAxe pipeline and script functions.
"""

import argparse
import os
from ast import literal_eval
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import numpy as np

from thoraxe import subexons
from thoraxe import utils
from thoraxe import transcript_info
from thoraxe.transcript_info import phases


def parse_command_line():
    """
    Parse command line.

    It uses argparse to parse thoraxe' command line arguments and returns
    the argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="thoraxe",
        description="""
        thoraxe is a tool to identify orthologous exonic regions (s-exons) in a
        set of transcripts from orthologous genes in a set of species.
        """,
        epilog="""
        It has been developed at LCQB (Laboratory of Computational and
        Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i',
        '--inputdir',
        help='Input directory. The input folder should have an Ensembl '
        'subfolder as the generated by transcript_query.',
        type=str,
        default='.')
    parser.add_argument(
        '-o',
        '--outputdir',
        help='Output directory, the indicated input directory is used '
        'by default.',
        type=str,
        default='')
    parser.add_argument('-a',
                        '--aligner',
                        help='Path to ProGraphMSA.',
                        type=str,
                        default='ProGraphMSA')
    parser.add_argument('-s',
                        '--maxtsl',
                        help='Maximum Transcript Support Level (TSL) to use '
                        'when TSL is available for a transcript.',
                        type=int,
                        default=3)
    parser.add_argument('-m',
                        '--minlen',
                        help='Minimum exon length.',
                        type=int,
                        default=4)
    parser.add_argument('-g',
                        '--mingenes',
                        help='Minimum number of genes to consider a path in '
                        'the splice graph.',
                        type=int,
                        default=1)
    parser.add_argument('-t',
                        '--mintranscripts',
                        help='Minimum number of transcripts to consider a '
                        'path in the splice graph.',
                        type=int,
                        default=2)
    parser.add_argument(
        '-c',
        '--coverage',
        help='Minimum alignment coverage of the shorter exon to include both '
        'exons in the same cluster.',
        type=float,
        default=80.0)
    parser.add_argument(
        '-p',
        '--identity',
        help='Minimum percent identity to include exons in the same cluster.',
        type=float,
        default=30.0)
    parser.add_argument('--gapopen',
                        help='Penalty for a gap opening.',
                        type=int,
                        default=-10)
    parser.add_argument('--gapextend',
                        help='Penalty for gap extensions.',
                        type=int,
                        default=-1)
    parser.add_argument('-r',
                        '--rescue_unaligned_subexons',
                        help='The sub-exons that do not align against any '
                        'other are deleted from their cluster, and they '
                        'could be reassigned during the sub-exon rescue '
                        'phase. By default, they are kept into their '
                        'original exon cluster.',
                        action='store_true')
    parser.add_argument(
        '--padding',
        help='Length of padding, Xs, in the chimeric alignment.',
        type=int,
        default=10)
    parser.add_argument(
        '-y',
        '--phylosofs',
        help='Save inputs to run PhyloSofS in the phylosofs folder.',
        action='store_true')
    parser.add_argument('--no_movements',
                        help='Do not move one/two residue sub-exon blocks.',
                        action='store_true')
    parser.add_argument('--no_disintegration',
                        help='Do not disintegrate one-residue-length s-exons.',
                        action='store_true')
    parser.add_argument(
        '--plot_chimerics',
        help='Save plotly/html plot for the chimeric alignments in the '
        '_intermediate folder.',
        action='store_true')
    parser.add_argument(
        '-l',
        '--specieslist',
        help='It could be a list of more than one species separated by commas '
        'and without spaces, e.g. homo_sapiens,mus_musculus, or a single file '
        'with the species list (one species per line). If nothing is '
        'indicated, all the available species are used.',
        default='')
    parser.add_argument(
        '--canonical_criteria',
        help='List of column names of the `path_table` separated by commas '
        'used to sort the row. If nothing is indicated, the following list is '
        'used: MinimumConservation,MinimumTranscriptWeightedConservation,'
        'MeanTranscriptWeightedConservation,TranscriptLength,TSL',
        default='')
    return parser


def get_transcripts(input_folder, max_tsl_level=3.0, species_list=None):
    """Return a DataFrame with the transcript information."""
    ensembl_info = os.path.join(input_folder, 'Ensembl')
    return transcript_info.read_transcript_info(
        os.path.join(ensembl_info, 'tsl.csv'),
        os.path.join(ensembl_info, 'exonstable.tsv'),
        os.path.join(ensembl_info, 'sequences.fasta'),
        max_tsl_level=max_tsl_level,
        remove_na=False,
        species_list=species_list)


def get_subexons(  # pylint: disable=too-many-arguments
        transcript_table,
        minimum_len=4,
        coverage_cutoff=80.0,
        percent_identity_cutoff=30.0,
        gap_open_penalty=10,
        gap_extend_penalty=1):
    """
    Return a DataFrame with subexon information and clustered exons.

    Exons are clustered according to their protein sequence using the blosum50
    matrix and the following arguments: minimum_len, coverage_cutoff,
    percent_identity_cutoff, gap_open_penalty and gap_extend_penalty.
    """
    clustered = transcript_info.exon_clustering(
        transcript_table,
        minimum_len=minimum_len,
        coverage_cutoff=coverage_cutoff,
        percent_identity_cutoff=percent_identity_cutoff,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty)
    subexon_table = subexons.create_subexon_table(clustered)
    return merge_clusters(subexon_table)


def _merge_clusters_in_list(cluster_lists):
    """Return a list of sets, each set has the clusters to merge."""
    output_list = []
    for input_list in cluster_lists:
        input_set = set(input_list)
        added = False
        for output_set in output_list:
            if output_set.intersection(input_set):
                output_set.update(input_set)
                added = True
        if not added:
            output_list.append(input_set)
    return output_list


def merge_clusters(subexon_table):
    """Merge 'Cluster's that share subexons."""
    clusters_df = subexon_table.groupby('SubexonIDCluster').agg(
        {'Cluster': lambda col: set(val for val in col if val != 0)})
    cluster_lists = clusters_df.loc[
        map(lambda x: len(x) > 1, clusters_df['Cluster']), 'Cluster'].apply(
            repr).drop_duplicates().apply(lambda x: sorted(literal_eval(x)))
    cluster_sets = _merge_clusters_in_list(cluster_lists)
    for cluster_set in cluster_sets:
        clusters = sorted(cluster_set)
        cluster_id = clusters[0]
        for cluster in clusters[1:]:
            cluster_index = subexon_table.index[subexon_table['Cluster'] ==
                                                cluster]
            for i in cluster_index:
                subexon_table.at[i, 'Cluster'] = cluster_id
    return subexon_table


def _outfile(output_folder, prefix, name, ext):
    """
    Return the output file name.

    >>> _outfile("", "gene_ids_", 10, ".txt")
    'gene_ids_10.txt'
    """
    return os.path.join(output_folder, prefix + str(name) + ext)


def _create_chimeric_msa(  # pylint: disable=too-many-arguments
        output_folder,
        cluster,
        subexon_df,
        gene2speciesname,
        connected_subexons,
        aligner='ProGraphMSA',
        padding='XXXXXXXXXX',
        species_list=None):
    """Return a modified subexon_df, the dict of chimerics and the msa."""
    subexon_df, subexon_matrix = subexons.alignment.create_subexon_matrix(
        subexon_df)
    chimerics = subexons.alignment.create_chimeric_sequences(
        subexon_df, subexon_matrix, connected_subexons, padding=padding)
    msa_file = _outfile(output_folder, "chimeric_alignment_", cluster,
                        ".fasta")
    if chimerics:
        chimerics = subexons.alignment.sort_species(chimerics,
                                                    gene2speciesname,
                                                    species_list)
        subexons.alignment.run_aligner(chimerics,
                                       aligner=aligner,
                                       output_path=msa_file)
        msa = subexons.alignment.read_msa_fasta(msa_file)
    else:
        if os.path.isfile(msa_file):
            os.remove(msa_file)
        msa = None
    return subexon_df, chimerics, msa


def _create_and_test_chimeric_msa(  # noqa pylint: disable=too-many-arguments,too-many-locals
        cluster2data,
        subexon_table,
        cluster_index,
        output_folder,
        cluster,
        gene2speciesname,
        connected_subexons,
        cutoff=30.0,
        min_col_number=4,
        aligner='ProGraphMSA',
        padding='XXXXXXXXXX',
        species_list=None,
        keep_single_subexons=False):
    """
    Update cluster2data and return a set of exons to delete.

    For each cluster, it adds the subexon dataframe, the chimeric sequences
    and the msa to cluster2data.
    """
    subexon_df = subexon_table.loc[cluster_index, :]
    subexon_df, chimerics, msa = _create_chimeric_msa(
        output_folder,
        cluster,
        subexon_df,
        gene2speciesname,
        connected_subexons,
        aligner=aligner,
        padding=padding,
        species_list=species_list)
    cluster2data[cluster] = (subexon_df, chimerics, msa)
    if msa is None:
        return {}
    return subexons.alignment.delete_subexons(
        subexon_df,
        chimerics,
        msa,
        cutoff=cutoff,
        min_col_number=min_col_number,
        keep_single_subexons=keep_single_subexons)


def create_chimeric_msa(  # pylint: disable=too-many-arguments,too-many-locals
        output_folder,
        subexon_table,
        gene2speciesname,
        connected_subexons,
        clusters=None,
        cutoff=30.0,
        min_col_number=4,
        aligner='ProGraphMSA',
        padding='XXXXXXXXXX',
        species_list=None,
        keep_single_subexons=False):
    """
    Return a dict from cluster to cluster data.

    For each cluster, there is a tuple with the subexon dataframe, the
    chimeric sequences and the msa.

    This function can take a clusters argument with the list of 'Cluster'
    identifiers to use. If that list is not given, all the positive 'Cluster'
    identifiers from subexon_table are used.
    """
    if clusters is None:
        clusters = {
            cluster
            for cluster in subexon_table['Cluster'] if cluster > 0
        }

    cluster2data = {}
    for cluster in clusters:
        cluster_index = subexon_table.index[subexon_table['Cluster'] ==
                                            cluster]
        to_delete = _create_and_test_chimeric_msa(
            cluster2data,
            subexon_table,
            cluster_index,
            output_folder,
            cluster,
            gene2speciesname,
            connected_subexons,
            cutoff=cutoff,
            min_col_number=min_col_number,
            aligner=aligner,
            padding=padding,
            species_list=species_list,
            keep_single_subexons=keep_single_subexons)
        while to_delete:
            delete = [
                subexon_table.loc[index, 'SubexonIDCluster'] in to_delete
                for index in cluster_index
            ]
            for index in cluster_index[delete]:
                subexon_table.at[index, 'Cluster'] = -1 * abs(
                    subexon_table.loc[index, 'Cluster'])
            cluster_index = cluster_index[np.logical_not(delete)]
            to_delete = _create_and_test_chimeric_msa(
                cluster2data,
                subexon_table,
                cluster_index,
                output_folder,
                cluster,
                gene2speciesname,
                connected_subexons,
                cutoff=cutoff,
                min_col_number=min_col_number,
                aligner=aligner,
                padding=padding,
                species_list=species_list,
                keep_single_subexons=keep_single_subexons)
    return cluster2data


def get_s_exons(  # noqa pylint: disable=too-many-arguments,too-many-locals
        output_folder,
        subexon_table,
        gene2speciesname,
        connected_subexons,
        minimum_len=4,
        coverage_cutoff=80.0,
        percent_identity_cutoff=30.0,
        gap_open_penalty=10,
        gap_extend_penalty=1,
        aligner='ProGraphMSA',
        padding='XXXXXXXXXX',
        movements=True,
        disintegration=True,
        species_list=None,
        keep_single_subexons=False):
    """Perform almost all the ThorAxe pipeline."""

    intermediate_output_path = utils.folders.create_subfolder(
        output_folder, '_intermediate')
    msa_output_path = utils.folders.create_subfolder(output_folder, 'msa')

    cluster2updated_data = {}
    cluster2data = create_chimeric_msa(
        intermediate_output_path,
        subexon_table,
        gene2speciesname,
        connected_subexons,
        cutoff=percent_identity_cutoff,
        min_col_number=minimum_len,
        aligner=aligner,
        padding=padding,
        species_list=species_list,
        keep_single_subexons=keep_single_subexons)
    modified_clusters = subexons.rescue.subexon_rescue_phase(
        cluster2data,
        subexon_table,
        minimum_len=minimum_len,
        coverage_cutoff=coverage_cutoff,
        percent_identity_cutoff=percent_identity_cutoff,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty,
        substitution_matrix=None)
    cluster2data.update(
        create_chimeric_msa(intermediate_output_path,
                            subexon_table,
                            gene2speciesname,
                            connected_subexons,
                            clusters=modified_clusters,
                            cutoff=percent_identity_cutoff,
                            min_col_number=minimum_len,
                            aligner=aligner,
                            padding=padding,
                            species_list=species_list,
                            keep_single_subexons=keep_single_subexons))

    for (cluster, (subexon_df, chimerics, msa)) in cluster2data.items():
        if (msa is not None) and chimerics:
            gene_ids = subexons.alignment.get_gene_ids(msa)
            msa_matrix = subexons.alignment.create_msa_matrix(chimerics, msa)

            if movements:
                (msa_numpy, msa_matrix
                 ) = subexons.alignment.move_problematic_block_clusters(
                     msa, msa_matrix)
            else:
                msa_numpy = np.array([list(rec) for rec in msa], np.character)

            if disintegration:
                msa_numpy, msa_matrix = subexons.alignment.disintegration(
                    msa_numpy, msa_matrix)

            if movements or disintegration:
                msa = MultipleSeqAlignment([
                    SeqRecord(Seq(row.tostring().decode('utf-8'),
                                  generic_protein),
                              id=gene_id)
                    for (gene_id, row) in zip(gene_ids, msa_numpy)
                ])
                # cluster2data[cluster] = (subexon_df, chimerics, msa)
            with open(
                    _outfile(intermediate_output_path, "gene_ids_", cluster,
                             ".txt"), 'w') as outfile:
                for item in gene_ids:
                    outfile.write("%s\n" % item)

            np.savetxt(_outfile(intermediate_output_path, "msa_matrix_",
                                cluster, ".txt"),
                       msa_matrix,
                       delimiter=",")

            colclusters = subexons.alignment.column_clusters(
                subexons.alignment.column_patterns(msa_matrix))

            sequences = subexons.alignment.msa2sequences(
                msa, gene_ids, padding)
            subexon_df = subexons.alignment.save_s_exons(
                subexon_df, sequences, gene_ids, colclusters, msa_output_path)
        else:
            msa = None
            gene_ids = None
            msa_matrix = None
            colclusters = None

        cluster2updated_data[cluster] = (subexon_df, chimerics, msa, gene_ids,
                                         msa_matrix, colclusters)
        subexon_df.to_csv(
            _outfile(intermediate_output_path, "subexon_table_", cluster,
                     ".csv"))

    return cluster2updated_data


def update_subexon_table(subexon_table, cluster2data):
    """Update the subexon table by adding the s-exon information."""
    columns_to_add = [
        'S_exon_Lengths', 'S_exon_Sequences', 'S_exons', 'SubexonIndex'
    ]
    key_col = 'SubexonIDCluster'
    subexon2data = {}
    for (_, data) in cluster2data.items():
        cluster_df = data[0]
        if cluster_df.shape[0] > 0:
            subexon2data.update(
                cluster_df.loc[:, [key_col] +
                               columns_to_add].drop_duplicates().set_index(
                                   key_col).to_dict('index'))

    columns = {colname: [] for colname in columns_to_add}
    for subexon in subexon_table[key_col]:
        for colname in columns_to_add:
            if subexon in subexon2data:
                columns[colname].append(subexon2data[subexon][colname])
            else:
                columns[colname].append('')

    for (colname, values) in columns.items():
        subexon_table[colname] = values

    return subexon_table


def add_s_exon_phases_and_coordinates(tbl):
    """
    Add s-exon genomic coordinates and phases to the tidy s-exon table.

    For each s-exon, it add the ``S_exon_Genomic_Sequence`` and their genomic
    coordinates, starting at ``1`` and defining the closed interval:
    ``[S_exon_CodingStart, S_exon_CodingEnd]``. Where ``S_exon_CodingStart``
    is greater or equal to ``S_exon_CodingEnd`` for genes in the reverse
    ``Strand``. It add also the ``S_exon_StartPhase`` and ``S_exon_EndPhase``
    for those intervals.

    An s-exon, defined at the protein level, can be an entire sub-exon or a
    part of it. In cases where an s-exon is identical to a sub-exon, the
    genomic coordinates and the phases are the same. However, when a sub-exon
    is split into several s-exons, coordinates and phases need to be calculated:

    **Example 1: s-exons in the positive strand.**

    The following sub-exon shares codons with the previous and following
    sub-exon. In this example, ``-`` represents introns and lower case
    characters represent exon nucleotides. The start and end phases are ``1``
    and ``2`` respectively:

    ``...v----vvwwwxxxyyyzz-----z...``

    Because we are in the positive strand, the residue ``V`` coming from the
    ``vvv`` shared codon is assigned to this sub-exon protein sequence. In the
    same way, the ``Z`` residue coming from the ``zzz`` shared exon is assigned
    to the beginning of the next sub-exon:

    ``VWXY``

    Let's say that this sub-exon is composed of two s-exons: ``VW`` and ``XY``.
    Then the corresponding genomic sequences are ``vvwww`` and ``xxxyyyzz`` and
    the ``(start, end)`` phases are ``(1, 0)`` and ``(0, 2)``.

    *Genomic coordinates.*

    Let's say that the first ``v`` has the genomic coordinate ``1``:

    ``v ---- vvwww xxxyyyzz ----- z``
    ``           1 11111111 12222 2``
    ``1 2345 67890 12345678 90123 4``

    In this way, the first s-exon has the genomic coordinates ``[6, 10]`` and
    the second has the genomic coordinates ``[11, 18]``. The first s-exon has
    ``2`` amino-acid residues, but ``5`` instead of ``6`` bases because the
    first codon is shared between two s-exons. In particular, the
    amino-acid residue ``V`` needs two bases in this s-exon to complete the
    codon started in the previous s-exon.

    **Example 2: s-exons in the negative strand.**

    ``v ---- vvwww xxxyyyzz ----- z``
    ``2 2222 11111 11111           ``
    ``4 3210 98765 43210987 65432 1``

    If the same sequence belongs to the negative strand, the sub-exon protein
    sequence is going to contain ``Z`` but not ``V``: ``WXYZ``. Therefore,
    possible s-exons are ``W`` and ``XYZ`` and the corresponding genomic
    sequences will be ``vvwww`` and ``xxxyyyzz`` and the ``(start, end)``
    phases will be ``(1, 0)`` and ``(0, 2)``.

    In this way, the first s-exon has the genomic coordinates ``[19, 15]`` and
    the second has the genomic coordinates ``[14, 7]``. The first s-exon has `
    `1`` amino-acid residue, but ``5`` instead of ``3`` bases because the first
    codon is shared between two s-exons. In particular, the amino-acid residue
    ``V`` needs two bases in this s-exon to complete the codon started in the
    previous s-exon.

    """
    tbl['S_exon_CodingStart'] = None
    tbl['S_exon_CodingEnd'] = None
    tbl['S_exon_StartPhase'] = None
    tbl['S_exon_EndPhase'] = None
    tbl['S_exon_Genomic_Sequence'] = None

    for _, group in tbl.groupby(['SubexonIDCluster', 'TranscriptIDCluster']):
        group_indices = group.index
        n_rows = len(group_indices)
        if n_rows > 1:
            end_coord = 0
            remaining_seq = ''
            for row_number, i in enumerate(group_indices):
                if row_number == 0:  # first s-exon
                    tbl.at[i,
                           'S_exon_StartPhase'] = tbl.loc[i,
                                                          'SubexonStartPhase']
                    tbl.at[i, 'S_exon_EndPhase'] = 0
                    (end_coord, s_exon_seq,
                     remaining_seq) = _get_first_s_exon(tbl.loc[i, :])
                    tbl.at[i, 'S_exon_CodingStart'] = tbl.loc[
                        i, 'SubexonCodingStart']
                    tbl.at[i, 'S_exon_CodingEnd'] = end_coord
                    tbl.at[i, 'S_exon_Genomic_Sequence'] = s_exon_seq
                elif row_number == (n_rows - 1):  # last s-exon
                    tbl.at[i, 'S_exon_StartPhase'] = 0
                    tbl.at[i, 'S_exon_EndPhase'] = tbl.loc[i,
                                                           'SubexonEndPhase']
                    tbl.at[i, 'S_exon_CodingStart'] = end_coord + \
                        tbl.loc[i, 'Strand']
                    tbl.at[i, 'S_exon_CodingEnd'] = tbl.loc[i,
                                                            'SubexonCodingEnd']
                    tbl.at[i, 'S_exon_Genomic_Sequence'] = remaining_seq

                else:
                    tbl.at[i, 'S_exon_StartPhase'] = 0
                    tbl.at[i, 'S_exon_EndPhase'] = 0
                    (start_coord, end_coord, s_exon_seq,
                     remaining_seq) = _get_internal_s_exon(
                         tbl.loc[i, :], end_coord, remaining_seq)
                    tbl.at[i, 'S_exon_CodingStart'] = start_coord
                    tbl.at[i, 'S_exon_CodingEnd'] = end_coord
                    tbl.at[i, 'S_exon_Genomic_Sequence'] = s_exon_seq
        elif n_rows == 1:
            i = group_indices[0]
            tbl.at[i, 'S_exon_CodingStart'] = tbl.loc[i, 'SubexonCodingStart']
            tbl.at[i, 'S_exon_CodingEnd'] = tbl.loc[i, 'SubexonCodingEnd']
            tbl.at[i, 'S_exon_StartPhase'] = tbl.loc[i, 'SubexonStartPhase']
            tbl.at[i, 'S_exon_EndPhase'] = tbl.loc[i, 'SubexonEndPhase']
            tbl.at[i, 'S_exon_Genomic_Sequence'] = tbl.loc[i,
                                                           'SubexonSequence']

    return tbl


def _get_first_s_exon(row):
    """
    Calculate the genomic end coordinate of the first s-exon

    Also returns the genomic nucleotide sequences of the s-exon and the rest of
    the sub-exon.

    >>> row = {'SubexonSequence': 'TTGCCGTCATGAGCAG'}
    >>> row['SubexonCodingStart'] = 2379656
    >>> row['SubexonStartPhase'] = 1
    >>> row['S_exon_Sequence'] = 'AV'
    >>> row['Strand'] = -1
    >>> _get_first_s_exon(row)
    (2379649, 'TTGCCGTC', 'ATGAGCAG')
    >>> row = {'SubexonSequence': 'CACAGGTGCAGCAGTGA'}
    >>> row['SubexonCodingStart'] = 48434879
    >>> row['SubexonStartPhase'] = 1
    >>> row['S_exon_Sequence'] = 'AQ'
    >>> row['Strand'] = 1
    >>> _get_first_s_exon(row)
    (48434883, 'CACAG', 'GTGCAGCAGTGA')
    >>> # SNAP25 ENSE00000906745
    >>> seq = 'AGTAAAGATGCTGGTATCAGGACTTTGGTTATGTTGGATGAACAAGGAG'
    >>> row = {'SubexonSequence': seq}
    >>> row['SubexonCodingStart'] = 10284724
    >>> row['SubexonStartPhase'] = 0
    >>> row['S_exon_Sequence'] = 'S'
    >>> row['Strand'] = 1
    >>> _get_first_s_exon(row)
    (10284726, 'AGT', 'AAAGATGCTGGTATCAGGACTTTGGTTATGTTGGATGAACAAGGAG')
    """
    genomic_seq = row['SubexonSequence']
    start_coordinates = row['SubexonCodingStart']
    start_phase = row['SubexonStartPhase']
    s_exon_len = len(row['S_exon_Sequence'])

    needed_bases = s_exon_len * 3

    bases_to_complete_previous_codon = phases.bases_to_complete_previous_codon(
        start_phase)
    if bases_to_complete_previous_codon:
        bases_in_previous_codon = 3 - bases_to_complete_previous_codon
    else:
        bases_in_previous_codon = 0

    if row['Strand'] == 1:
        total_len = needed_bases - bases_in_previous_codon
        end_coordinates = start_coordinates + total_len - 1
    else:
        total_len = needed_bases + bases_to_complete_previous_codon
        end_coordinates = start_coordinates - total_len + 1

    return end_coordinates, genomic_seq[:total_len], genomic_seq[total_len:]


def _get_internal_s_exon(row, previous_end, seq):
    """
    Calculate the genomic coordinates of internal s-exons

    Also returns the genomic nucleotide sequences of the s-exon and the rest of
    the sub-exon.

    >>> row = {'S_exon_Sequence': 'MS'}
    >>> row['Strand'] = -1
    >>> _get_internal_s_exon(row, 2379649, 'ATGAGCAG')
    (2379648, 2379643, 'ATGAGC', 'AG')
    >>> row['Strand'] = 1
    >>> _get_internal_s_exon(row, 2379649, 'ATGAGCAG')
    (2379650, 2379655, 'ATGAGC', 'AG')
    """
    total_len = len(row['S_exon_Sequence']) * 3
    if row['Strand'] == 1:
        start_coordinates = previous_end + 1
        end_coordinates = previous_end + total_len
    else:
        start_coordinates = previous_end - 1
        end_coordinates = previous_end - total_len
    return start_coordinates, end_coordinates, seq[:total_len], seq[total_len:]


def get_s_exon_msas(output_folder):
    """
    Return a dict of the s_exon MSAs.
    """
    s_exon_msas = dict()
    msa_path = os.path.join(output_folder, 'msa')
    msa_files = os.listdir(msa_path)
    for msa_file in msa_files:
        if 'msa_s_exon_' in msa_file and '.fasta' in msa_file:
            fields = msa_file.split('_')
            s_exon = '{}_{}'.format(fields[3], fields[4].split('.')[0])
            s_exon_msas[s_exon] = subexons.alignment.read_msa_fasta(
                os.path.join(msa_path, msa_file))
    return s_exon_msas


def main():  # pylint: disable=too-many-locals
    """Perform Pipeline."""
    args = parse_command_line().parse_args()

    species_list = utils.species.get_species_list(args.specieslist)

    input_folder = os.path.abspath(args.inputdir)
    output_folder = input_folder if args.outputdir == '' else args.outputdir
    output_folder = utils.folders.create_subfolder(output_folder, 'thoraxe')

    intermediate_output_path = utils.folders.create_subfolder(
        output_folder, '_intermediate')

    transcript_table = get_transcripts(input_folder,
                                       float(args.maxtsl),
                                       species_list=species_list)
    subexon_table = get_subexons(transcript_table,
                                 minimum_len=args.minlen,
                                 coverage_cutoff=args.coverage,
                                 percent_identity_cutoff=args.identity,
                                 gap_open_penalty=args.gapopen,
                                 gap_extend_penalty=args.gapextend)

    transcript_table.to_csv(
        os.path.join(intermediate_output_path, "transcript_table.csv"))

    gene2speciesname = subexons.alignment.gene2species(transcript_table)
    connected_subexons = subexons.alignment.subexon_connectivity(subexon_table)

    cluster2data = get_s_exons(
        output_folder,
        subexon_table,
        gene2speciesname,
        connected_subexons,
        minimum_len=args.minlen,
        coverage_cutoff=args.coverage,
        percent_identity_cutoff=args.identity,
        gap_open_penalty=args.gapopen,
        gap_extend_penalty=args.gapextend,
        aligner=args.aligner,
        padding='X' * args.padding,
        movements=not args.no_movements,
        disintegration=not args.no_disintegration,
        species_list=species_list,
        keep_single_subexons=not args.rescue_unaligned_subexons)

    s_exon_msas = get_s_exon_msas(output_folder)

    subexon_table.to_csv(
        os.path.join(intermediate_output_path, "subexon_table.csv"))

    if args.plot_chimerics:
        subexons.plot.plot_msa_subexons(cluster2data, intermediate_output_path)

    subexon_table = update_subexon_table(subexon_table, cluster2data)
    subexon_table = subexons.alignment.impute_missing_s_exon(subexon_table)

    tidy_table = subexons.tidy.get_tidy_table(subexon_table, gene2speciesname)
    add_s_exon_phases_and_coordinates(tidy_table)
    tidy_table.to_csv(os.path.join(output_folder, "s_exon_table.csv"),
                      index=False)

    if args.phylosofs:
        s_exon_2_char = subexons.phylosofs.phylosofs_inputs(
            tidy_table, os.path.join(input_folder, 'Ensembl'), output_folder)
    else:
        s_exon_2_char = {}

    (node2genes, edge2genes, node2transcripts, edge2transcripts, edge2trx_cons
     ) = subexons.graph.nodes_and_edges2genes_and_transcripts(tidy_table)

    splice_graph_filename = os.path.join(output_folder, "splice_graph.gml")

    subexons.graph.splice_graph_gml(splice_graph_filename, s_exon_msas,
                                    node2genes, edge2genes, node2transcripts,
                                    edge2transcripts, edge2trx_cons,
                                    s_exon_2_char)

    path_table, ases_table = subexons.ases.conserved_ases(
        tidy_table,
        transcript_table,
        s_exon_msas,
        splice_graph_filename,
        min_genes=args.mingenes,
        min_transcripts=args.mintranscripts,
        column_order=args.canonical_criteria.split(",")
        if args.canonical_criteria else None)

    path_table.to_csv(os.path.join(output_folder, "path_table.csv"),
                      index=False)
    ases_table.to_csv(os.path.join(output_folder, "ases_table.csv"),
                      index=False)


if __name__ == '__main___':
    main()
