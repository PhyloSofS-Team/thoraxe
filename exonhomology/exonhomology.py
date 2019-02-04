import argparse
import logging
import os
import re
import numpy as np
from . import transcript_info
from . import subexons


def parse_command_line():
    """
    Parse command line.

    It uses argparse to parse exonhomology' command line arguments and returns
    the argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="exonhomology",
        description="""
        exonhomology is a tool to identify homologous exons in a set of
        transcripts from the same gene and different species.
        """,
        epilog="""
        It has been developed at LCQB (Laboratory of Computational and
        Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
        """)
    parser.add_argument(
        '-i', '--inputdir', help='input directory', type=str, default='.')
    parser.add_argument(
        '-g',
        '--genename',
        help='gene name using the transcript_query output format '
        '(e.g. MAPK8_ENSG00000107643). By default, the name of the last '
        'directory in --inputdir (-i) is used',
        type=str,
        default='')
    parser.add_argument(
        '-o', '--outputdir', help='output directory', type=str, default='.')
    parser.add_argument(
        '-t', '--mafft_path', help='path to MAFFT', type=str, default='mafft')
    parser.add_argument(
        '-m', '--minlen', help='minimum exon length', type=int, default=4)
    parser.add_argument(
        '-c',
        '--coverage',
        help='minimum alignment coverage of the shorter exon to include both '
        'exons in the same cluster',
        type=float,
        default=80.0)
    parser.add_argument(
        '-p',
        '--identity',
        help='minimum percent identity to include exons in the same cluster',
        type=float,
        default=30.0)
    parser.add_argument(
        '--gapopen', help='penalty for a gap opening', type=int, default=10)
    parser.add_argument(
        '--gapextend', help='penalty for gap extensions', type=int, default=1)

    return parser.parse_args()


def _get_gene_name(input_folder, gene_name=None):
    r"""
    Return the input_folder and the gene_name.

    It gives a warning if the gene_name hasn't the correct format.

    >>> _get_gene_name('exonhomology/tests/data/MAPK8_ENSG00000107643')
    ('exonhomology/tests/data/MAPK8_ENSG00000107643', 'MAPK8_ENSG00000107643')
    >>> _get_gene_name('exonhomology/tests/data/MAPK8')
    ('exonhomology/tests/data/MAPK8', 'MAPK8')
    """
    input_folder = os.path.abspath(input_folder)
    if gene_name is None:
        gene_name = os.path.basename(os.path.normpath(input_folder))
    # \.?[0-9]* is needed to accept version numbers
    if re.match(r'^\w+_ENS\w+$', gene_name) is None:
        logging.warning(
            'The gene_name %s does not have the proper format. %s', gene_name,
            'It should be the gene name followed by _ and the '
            'Ensembl stable ID (e.g. MAPK8_ENSG00000107643).')
    return input_folder, gene_name


def get_transcripts(input_folder, gene_name):
    """Return a DataFrame with the transcript information."""
    return transcript_info.read_transcript_info(
        os.path.join(input_folder, 'TSL', gene_name + '_TSL.csv'),
        os.path.join(input_folder, 'TablesExons',
                     gene_name + '_exonstable.tsv'),
        os.path.join(input_folder, 'Sequences', gene_name + '.fasta'),
        remove_na=False)


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
    return subexons.create_subexon_table(clustered)


def _outfile(output_folder, prefix, name, ext):
    """
    Return the output file name.

    >>> _outfile("", "gene_ids_", 10, ".txt")
    'gene_ids_10.txt'
    """
    return os.path.join(output_folder, prefix + str(name) + ext)


def get_homologous_subexons(  # noqa pylint: disable=too-many-arguments,too-many-locals
        outdir,
        name,
        subexon_df,
        gene2speciesname,
        connected_subexons,
        mafft_path='mafft'):
    """Perform almost all the pipeline."""
    subexon_df, subexon_matrix = subexons.alignment.create_subexon_matrix(
        subexon_df)
    chimerics = subexons.alignment.create_chimeric_sequences(
        subexon_df, subexon_matrix, connected_subexons)
    chimerics = subexons.alignment.sort_species(chimerics, gene2speciesname)

    msa_file = _outfile(outdir, "chimeric_alignment_", name, ".fasta")
    subexons.alignment.run_mafft(
        chimerics, mafft_path=mafft_path, output_path=msa_file)

    if os.path.getsize(msa_file) > 0:
        gene_ids, msa_matrix = subexons.alignment.create_msa_matrix(
            chimerics, msa_file)

        with open(_outfile(outdir, "gene_ids_", name, ".txt"), 'w') as outfile:
            for item in gene_ids:
                outfile.write("%s\n" % item)

        np.savetxt(
            _outfile(outdir, "msa_matrix_", name, ".txt"),
            msa_matrix,
            delimiter=",")

        colclusters = subexons.alignment.column_clusters(
            subexons.alignment.column_patterns(msa_matrix))

        if subexons.plot._PLOT:  # pylint: disable=protected-access
            subexons.plot.plot_msa_subexons(
                gene_ids,
                msa_matrix,
                outfile=_outfile(outdir, "chimeric_alignment_", name, ".png"),
                subexon_table=subexon_df)

    else:
        gene_ids = None
        msa_matrix = None
        colclusters = None

    subexon_df.to_csv(_outfile(outdir, "subexon_table_", name, ".csv"))

    return (subexon_df, chimerics, subexon_matrix, gene_ids, msa_matrix,
            msa_file, colclusters)


def main():
    """Perform Pipeline."""
    args = parse_command_line()

    input_folder, gene_name = _get_gene_name(
        args.inputdir,
        gene_name=args.genename if args.genename != '' else None)
    transcript_table = get_transcripts(input_folder, gene_name)
    subexon_table = get_subexons(
        transcript_table,
        minimum_len=args.minlen,
        coverage_cutoff=args.coverage,
        percent_identity_cutoff=args.identity,
        gap_open_penalty=args.gapopen,
        gap_extend_penalty=args.gapextend)

    transcript_table.to_csv(
        os.path.join(args.outputdir, "transcript_table.csv"))
    subexon_table.to_csv(os.path.join(args.outputdir, "subexon_table.csv"))

    gene2speciesname = subexons.alignment.gene2species(transcript_table)
    connected_subexons = subexons.alignment.subexon_connectivity(subexon_table)

    for name, subgroup in subexon_table.groupby('Cluster'):
        if subgroup.shape[0] > 1:
            get_homologous_subexons(
                args.outputdir,
                name,
                subgroup,
                gene2speciesname,
                connected_subexons,
                mafft_path=args.mafft_path)


if __name__ == '__main___':
    main()
