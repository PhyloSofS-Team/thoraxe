"""plot: Generate js data for plotting the MSA matrix using plotly."""

import codecs
import os
import json
from shutil import copyfile
import pkg_resources
import numpy as np
import pandas as pd

_HTML_FILE = pkg_resources.resource_filename(
    'thoraxe', 'subexons/assets/cluster_plots.html')


def _constitutive_value(values):
    """
    Return 0 for alternative, 1 for constitutive.

    It returns 0.33 if there is only one transcript, and 0.66 if the subexon
    is present in all the transcripts except 1.
    """
    n_trx_gene, n_trx_subexon, fraction = values
    if np.isnan(fraction):
        return np.nan
    if n_trx_gene == 1:
        return 0.33
    if fraction == 1.0:
        return 1.0
    if n_trx_gene - 1 == n_trx_subexon:
        return 0.66
    return 0.0


def create_python_structure(cluster2data):
    """Python structured data to be translated to JavaScript."""
    structured_data = {}
    for cluster, data in cluster2data.items():
        msa = data[2]
        if msa is None:
            continue

        subexon_table = data[0]
        gene_ids = data[3]
        msa_matrix = data[4]

        fraction_data = pd.Series(
            zip(subexon_table['TranscriptsInGene'].values,
                subexon_table['TranscriptsWithSubexon'].values,
                subexon_table['TranscriptFraction'].values),
            index=subexon_table['SubexonIndex']).to_dict()
        z_constitutive = [[
            _constitutive_value(
                fraction_data.get(subexon, (np.nan, np.nan, np.nan)))
            for subexon in row
        ] for row in msa_matrix]

        structured_data['cluster_{}'.format(cluster)] = {
            'x': list(range(0, msa.get_alignment_length())),
            'y': gene_ids,
            'z_cluster': msa_matrix.tolist(),
            'z_constitutive': z_constitutive,
            'text': [list(record) for record in msa]
        }

    return structured_data


def plot_msa_subexons(cluster2data, output_folder):
    """Save html and js data in the output_folder."""
    structured_data = create_python_structure(cluster2data)

    js_file = os.path.join(output_folder, 'cluster_data.js')
    with open(js_file, 'w') as file:
        file.write('var clusterData = ')

    json.dump(
        structured_data,
        codecs.open(js_file, 'a', encoding='utf-8'),
        separators=(',', ':'),
        sort_keys=True,
        indent=4)

    with open(js_file, 'a') as file:
        file.write(';')

    copyfile(_HTML_FILE, os.path.join(output_folder, 'cluster_plots.html'))
