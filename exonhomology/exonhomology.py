import logging
import os
import re
# import transcript_info
# import subexons


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

    #
    # trx_data = transcript_info.read_transcript_info(
    # os.path.join(input_folder, 'TSL', folder + '_TSL.csv'),
    # os.path.join(input_folder, 'TablesExons', folder + '_exonstable.tsv'),
    # os.path.join(input_folder, 'Sequences', folder + '.fasta'),
    # remove_na=False)
    # data = transcript_info.exon_clustering(trx_data)
    # return subexons.create_subexon_table(data), trx_data
    #
    #
    # # CAMK2A_ENSG00000070808 GRIN1_ENSG00000169258  MAPK8_ENSG00000107643
    # data, trx_data = _get_clustered_trx_data('MAPK8_ENSG00000107643')
