import logging
import os
import re
import numpy as np
from . import transcript_info
from . import subexons


def _get_gene_name(input_folder, gene_name=None):
    r"""
    Return the input_folder and the gene_name.

    It gives a warning if the gene_name hasn't the correct format.

    >>> _get_gene_name('exonhomology/tests/data/MAPK8_ENSG00000107643')
    ('exonhomology/tests/data/MAPK8_ENSG00000107643', 'MAPK8_ENSG00000107643')
    >>> _get_gene_name('exonhomology/tests/data/MAPK8')
    ('exonhomology/tests/data/MAPK8', 'MAPK8')
    """
    if gene_name is None:
        gene_name = os.path.basename(os.path.normpath(input_folder))
    # \.?[0-9]* is needed to accept version numbers
    if re.match(r'^\w+_ENS\w+$', gene_name) is None:
        logging.warning(
            'The gene_name %s does not have the proper format. %s', gene_name,
            'It should be the gene name followed by _ and the '
            'Ensembl stable ID (e.g. MAPK8_ENSG00000107643).')
    return input_folder, gene_name


def _get_transcripts(input_folder, gene_name):
    """Return a DataFrame with the transcript information."""
    return transcript_info.read_transcript_info(
        os.path.join(input_folder, 'TSL', gene_name + '_TSL.csv'),
        os.path.join(input_folder, 'TablesExons',
                     gene_name + '_exonstable.tsv'),
        os.path.join(input_folder, 'Sequences', gene_name + '.fasta'),
        remove_na=False)


def _get_subexons(  # pylint: disable=too-many-arguments
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


def _get_homologous_subexons(outdir, name, subexon_df, gene2speciesname,
                             connected_subexons):
    """Perform almost all the pipeline."""
    subexon_df, subexon_matrix = subexons.alignment.create_subexon_matrix(
        subexon_df)
    chimerics = subexons.alignment.create_chimeric_sequences(
        subexon_df, subexon_matrix, connected_subexons)
    chimerics = subexons.alignment.sort_species(chimerics, gene2speciesname)

    msa_file = _outfile(outdir, "chimeric_alignment_", name, ".fasta")
    subexons.alignment.run_mafft(chimerics, output_path=msa_file)
    gene_ids, msa_matrix = subexons.alignment.create_msa_matrix(
        chimerics, msa_file)

    with open(_outfile(outdir, "gene_ids_", name, ".txt"), 'w') as outfile:
        for item in gene_ids:
            outfile.write("%s\n" % item)

    subexons.plot.plot_msa_subexons(
        gene_ids,
        msa_matrix,
        outfile=_outfile(outdir, "chimeric_alignment_", name, ".png"),
        subexon_table=subexon_df)

    np.savetxt(
        _outfile(outdir, "msa_matrix_", name, ".txt"),
        msa_matrix,
        delimiter=",")

    regions = subexons.alignment.minimal_regions(msa_matrix)
    with open(_outfile(outdir, "msa_subexon_regions_", name, ".txt"),
              'w') as regfile:
        for region in regions:
            regfile.write('%s\n' % str(region))

    subexon_df.to_csv(_outfile(outdir, "subexon_table_", name, ".csv"))
