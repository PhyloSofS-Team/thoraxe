"""thoraxe pipeline functions."""

import argparse
import os
from ast import literal_eval
import numpy as np
from thoraxe import transcript_info
from thoraxe import subexons
from thoraxe import utils


def parse_command_line():
    """
    Parse command line.

    It uses argparse to parse thoraxe' command line arguments and returns
    the argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="thoraxe",
        description="""
        thoraxe is a tool to identify homologous exons in a set of
        transcripts from the same gene and different species.
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
                        help='Path to Clustal Omega.',
                        type=str,
                        default='clustalo')
    parser.add_argument('-m',
                        '--minlen',
                        help='Minimum exon length.',
                        type=int,
                        default=4)
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
                        default=10)
    parser.add_argument('--gapextend',
                        help='Penalty for gap extensions.',
                        type=int,
                        default=1)
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
    parser.add_argument(
        '--plot_chimerics',
        help='Save plotly/html plot for the chimeric alignments in the '
        '_intermediate_outputs folder.',
        action='store_true')
    parser.add_argument(
        '-l',
        '--specieslist',
        help='It could be a list of more than one species separated by commas '
        'and without spaces, e.g. homo_sapiens,mus_musculus, or a single file '
        'with the species list (one species per line). If nothing is '
        'indicated, all the available species are used.',
        default='')

    return parser.parse_args()


def get_transcripts(input_folder, species_list=None):
    """Return a DataFrame with the transcript information."""
    ensembl_info = os.path.join(input_folder, 'Ensembl')
    return transcript_info.read_transcript_info(
        os.path.join(ensembl_info, 'tsl.csv'),
        os.path.join(ensembl_info, 'exonstable.tsv'),
        os.path.join(ensembl_info, 'sequences.fasta'),
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


def merge_clusters(subexon_table):
    """Merge 'Cluster's that share subexons."""
    clusters_df = subexon_table.groupby('SubexonIDCluster').agg(
        {'Cluster': lambda col: set(val for val in col if val != 0)})
    cluster_lists = clusters_df.loc[map(lambda x: len(
        x) > 1, clusters_df['Cluster']), 'Cluster'].apply(
            repr).drop_duplicates().apply(lambda x: sorted(literal_eval(x)))
    for clusters in cluster_lists:
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
        aligner='clustalo',
        padding='XXXXXXXXXX',
        species_list=None):
    """Return a modified subexon_df, the dict of chimerics and the msa."""
    subexon_df, subexon_matrix = subexons.alignment.create_subexon_matrix(
        subexon_df)
    chimerics = subexons.alignment.create_chimeric_sequences(
        subexon_df, subexon_matrix, connected_subexons, padding=padding)
    chimerics = subexons.alignment.sort_species(chimerics, gene2speciesname,
                                                species_list)
    msa_file = _outfile(output_folder, "chimeric_alignment_", cluster,
                        ".fasta")
    subexons.alignment.run_aligner(chimerics,
                                   aligner=aligner,
                                   output_path=msa_file)
    msa = subexons.alignment.read_msa_fasta(msa_file)
    return subexon_df, chimerics, msa


def _create_and_test_chimeric_msa(  # pylint: disable=too-many-arguments
        cluster2data,
        subexon_table,
        cluster_index,
        output_folder,
        cluster,
        gene2speciesname,
        connected_subexons,
        cutoff=30.0,
        min_col_number=4,
        aligner='clustalo',
        padding='XXXXXXXXXX',
        species_list=None):
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
    return subexons.alignment.delete_subexons(subexon_df,
                                              chimerics,
                                              msa,
                                              cutoff=cutoff,
                                              min_col_number=min_col_number)


def create_chimeric_msa(  # pylint: disable=too-many-arguments,too-many-locals
        output_folder,
        subexon_table,
        gene2speciesname,
        connected_subexons,
        clusters=None,
        cutoff=30.0,
        min_col_number=4,
        aligner='clustalo',
        padding='XXXXXXXXXX',
        species_list=None):
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
            species_list=species_list)
        while to_delete:
            delete = [
                subexon_table.loc[index, 'SubexonIDCluster'] in to_delete
                for index in cluster_index
            ]
            for index in cluster_index[delete]:
                subexon_table.at[index, 'Cluster'] = -1 * subexon_table.loc[
                    index, 'Cluster']
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
                species_list=species_list)
    return cluster2data


def get_homologous_subexons(  # noqa pylint: disable=too-many-arguments,too-many-locals
        output_folder,
        subexon_table,
        gene2speciesname,
        connected_subexons,
        minimum_len=4,
        coverage_cutoff=80.0,
        percent_identity_cutoff=30.0,
        gap_open_penalty=10,
        gap_extend_penalty=1,
        aligner='clustalo',
        padding='XXXXXXXXXX',
        species_list=None):
    """Perform almost all the pipeline."""

    intermediate_output_path = utils.folders.create_subfolder(
        output_folder, '_intermediate_outputs')

    cluster2updated_data = {}
    cluster2data = create_chimeric_msa(intermediate_output_path,
                                       subexon_table,
                                       gene2speciesname,
                                       connected_subexons,
                                       cutoff=percent_identity_cutoff,
                                       min_col_number=minimum_len,
                                       aligner=aligner,
                                       padding=padding,
                                       species_list=species_list)
    modified_clusters = subexons.rescue.subexon_rescue_phase(
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
                            species_list=species_list))

    for (cluster, (subexon_df, chimerics, msa)) in cluster2data.items():
        if msa is not None:
            gene_ids = subexons.alignment.get_gene_ids(msa)
            msa_matrix = subexons.alignment.create_msa_matrix(chimerics, msa)

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
            subexon_df = subexons.alignment.save_homologous_subexons(
                subexon_df, sequences, gene_ids, colclusters, output_folder)
        else:
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
    """Update the subexon table by adding the homologous exon information."""
    columns_to_add = [
        'HomologousExonLengths', 'HomologousExonSequences', 'HomologousExons',
        'SubexonIndex'
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


def main():
    """Perform Pipeline."""
    args = parse_command_line()

    species_list = utils.species.get_species_list(args.specieslist)

    input_folder = os.path.abspath(args.inputdir)
    output_folder = input_folder if args.outputdir == '' else args.outputdir

    intermediate_output_path = utils.folders.create_subfolder(
        output_folder, '_intermediate_outputs')

    transcript_table = get_transcripts(input_folder, species_list=species_list)
    subexon_table = get_subexons(transcript_table,
                                 minimum_len=args.minlen,
                                 coverage_cutoff=args.coverage,
                                 percent_identity_cutoff=args.identity,
                                 gap_open_penalty=args.gapopen,
                                 gap_extend_penalty=args.gapextend)

    transcript_table.to_csv(
        os.path.join(intermediate_output_path, "transcript_table.csv"))
    subexon_table.to_csv(
        os.path.join(intermediate_output_path, "subexon_table.csv"))

    gene2speciesname = subexons.alignment.gene2species(transcript_table)
    connected_subexons = subexons.alignment.subexon_connectivity(subexon_table)

    cluster2data = get_homologous_subexons(
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
        species_list=species_list)

    if args.plot_chimerics:
        subexons.plot.plot_msa_subexons(cluster2data, intermediate_output_path)

    subexon_table = update_subexon_table(subexon_table, cluster2data)
    subexon_table = subexons.alignment.impute_missing_orthologous_exon_group(
        subexon_table)
    subexon_table.to_csv(os.path.join(output_folder,
                                      "homologous_subexons.csv"))

    node2genes, edge2genes = subexons.graph.nodes_and_edges2genes(
        subexon_table)
    subexons.graph.splice_graph_gml(
        os.path.join(output_folder, "splice_graph.gml"), node2genes,
        edge2genes)

    if args.phylosofs:
        subexons.phylosofs.phylosofs_inputs(
            subexon_table, os.path.join(input_folder, 'Ensembl'),
            output_folder)

    tidy_table = subexons.tidy.get_tidy_table(subexon_table, gene2speciesname)
    tidy_table.to_csv(os.path.join(output_folder, "homologous_exon_table.csv"))


if __name__ == '__main___':
    main()
