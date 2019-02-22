"""plot: Plot the MSA matrix using plotly."""

import numpy as np
import pandas as pd
import plotly.figure_factory as ff
from plotly.offline import plot


def _create_heatmap(color_values, text_values, gene_ids):
    """Create a good looking and interactive plotly heatmap."""
    figure = ff.create_annotated_heatmap(
        color_values, y=gene_ids, annotation_text=text_values, xgap=2, ygap=2)
    layout = figure['layout']
    xaxis = layout['xaxis']
    yaxis = layout['yaxis']
    xaxis['showgrid'] = False
    xaxis['zeroline'] = False
    yaxis['showgrid'] = False
    yaxis['zeroline'] = False
    yaxis['automargin'] = True
    yaxis['autorange'] = 'reversed'
    return figure


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


def plot_msa_subexons(msa,
                      msa_matrix,
                      subexon_table=None,
                      outfile='alignment.html'):
    """
    Save a plot from the msa matrix, color indicates subexons.

    If `subexon_table` is `None`, the binary plot with the constitutive
    exons is not generated. Otherwise, it'll be in the `_constitutive.html`
    file.
    """
    assert outfile.endswith('.html')
    seq_ids = [seq.id for seq in msa]
    figure = _create_heatmap(msa_matrix, msa, seq_ids)
    figure['data'][0]['colorscale'] = 'Portland'
    plot(figure, filename=outfile, auto_open=False)
    if subexon_table is not None:
        fraction_data = pd.Series(
            zip(subexon_table['Transcripts in gene'].values,
                subexon_table['Number of transcripts for subexon'].values,
                subexon_table['Transcript fraction'].values),
            index=subexon_table['SubexonIndex']).to_dict()

        figure['data'][0]['z'] = [[
            _constitutive_value(
                fraction_data.get(subexon, (np.nan, np.nan, np.nan)))
            for subexon in row
        ] for row in msa_matrix]

        figure['data'][0]['colorscale'] = [
            [0.2, 'rgb(167, 219, 148)'],  # light
            [0.4, 'rgb(76, 147, 50)'],  # 0.33
            [0.7, 'rgb(44, 111, 19)'],  # 0.66
            [1.0, 'rgb(21, 74, 1)']  # dark
        ]
        plot(
            figure,
            filename=outfile.replace('.html', '_constitutive.html'),
            auto_open=False)
