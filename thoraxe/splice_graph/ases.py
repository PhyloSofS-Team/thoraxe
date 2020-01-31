"""
ases: Function to detect conserved alternative splicing events.

It finds the canonical path in the splice graph to detect conserved ASEs.
"""
import networkx as nx
import pandas as pd


def _get_seq(string):
    """
    Return the input string or an empty string with np.nan.

    >>> _get_seq('QVQQ')
    'QVQQ'
    >>> import numpy as np
    >>> _get_seq(np.nan)
    ''
    """
    if isinstance(string, str):
        return string
    return ""


def _get_transcript_scores(table, graph):
    """
    Return a DataFrame with the needed data to choose the canonical path.
    """
    data = {
        'TranscriptIDCluster': [],
        'TranscriptLength': [],
        'MinimumTranscriptWeightedConservation': [],
        'IsHuman': [],
        'Path': []
    }
    for (trx, subdf) in table.groupby('TranscriptIDCluster'):
        n_rows = len(subdf)
        s_exon_len = [len(_get_seq(subdf.S_exon_Sequence.iloc[0]))]
        path = ["start", subdf.S_exonID.iloc[0]]
        score = [1.0]  # 1.0 is the maximum and we select the minimum
        if n_rows >= 2:
            for i in range(1, n_rows):
                s_exon_1 = subdf.S_exonID.iloc[i - 1]
                s_exon_2 = subdf.S_exonID.iloc[i]
                s_exon_len.append(len(_get_seq(subdf.S_exon_Sequence.iloc[i])))
                path.append(subdf.S_exonID.iloc[i])
                score.append(
                    graph.get_edge_data(
                        s_exon_1,
                        s_exon_2)['transcript_weighted_conservation'])
        path.append("stop")
        data['TranscriptIDCluster'].append(trx)
        data['TranscriptLength'].append(sum(s_exon_len))
        data['MinimumTranscriptWeightedConservation'].append(min(score))
        data['IsHuman'].append(int(subdf.Species.iloc[0] == 'homo_sapiens'))
        data['Path'].append(",".join(path))
    data_frame = pd.DataFrame(data)
    data_frame.sort_values([
        'MinimumTranscriptWeightedConservation', 'IsHuman', 'TranscriptLength'
    ],
                           ascending=False,
                           inplace=True)
    return data_frame


def get_canonical_path(table, graph):
    """
    Return the canonical path.
    """
    trx_data = _get_transcript_scores(table, graph)
    path = trx_data.Path.iloc[0]
    return path.split(',')
