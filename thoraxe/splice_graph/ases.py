"""
ases: Function to detect conserved alternative splicing events.

It finds the canonical path in the splice graph to detect conserved ASEs.
"""
import collections
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


def _n_genes(graph, src, dst):
    """
    Return the number of genes showing a particular edge (src, dst).
    """
    return graph.get_edge_data(src, dst)['genes'].count(',') + 1


def _get_sources_and_destinies(graph, canonical_path, min_genes=1):
    """
    Return a list of sources and destinies in the canonical path.
    """
    sources = collections.OrderedDict()
    destinies = collections.OrderedDict()
    for (i, node) in enumerate(canonical_path):

        outneighbors = set(graph.successors(node))
        if node != 'stop':
            outneighbors.remove(canonical_path[i + 1])
        if sum(
                _n_genes(graph, node, dst) >= min_genes
                for dst in outneighbors) >= 1:
            sources[node] = i

        in_neighbors = set(graph.predecessors(node))
        if node != 'start':
            in_neighbors.remove(canonical_path[i - 1])
        if sum(
                _n_genes(graph, src, node) >= min_genes
                for src in in_neighbors) >= 1:
            destinies[node] = i
    return sources, destinies


def _is_alternative_path(  # pylint: disable=too-many-arguments
        graph,
        path,
        canonical_set,
        i,
        j,
        min_genes=1):
    """
    True if path is a good alternative path to the canonical path in the graph.
    """
    internal_nodes = path[1:-1]

    if (not internal_nodes) and (j - i == 1):
        return False

    for node in internal_nodes:
        if node in canonical_set:
            return False

    n_nodes = len(path)
    if n_nodes == 2:
        return _n_genes(graph, path[0], path[1]) >= min_genes

    genes = set(graph.get_edge_data(path[0], path[1])['genes'].split(','))
    if len(genes) >= min_genes:
        for k in range(2, n_nodes):
            genes.intersection_update(
                graph.get_edge_data(path[k - 1], path[k])['genes'].split(','))
            if len(genes) < min_genes:
                return False
    else:
        return False

    return len(genes) >= min_genes


def _is_mutually_exclusive(graph, canonical_internal, alternative_internal):
    """
    Return true if the is no edges between the paths.
    """
    for i in canonical_internal:
        for j in alternative_internal:
            if graph.has_edge(i, j) or graph.has_edge(j, i):
                return False
    return True


def _define_aes(graph, canonical_path, i, j, alternative_path):
    """
    Return a string indicating the kind of alternative splicing event.
    """
    canonical_internal = canonical_path[i + 1:j]
    alternative_internal = alternative_path[1:-1]
    if canonical_path[i] == 'start':
        return "alternative_start"
    if canonical_path[j] == 'stop':
        return "alternative_end"
    if not canonical_internal:
        return "insertion"
    if not alternative_internal:
        return "deletion"
    if _is_mutually_exclusive(graph, canonical_internal, alternative_internal):
        return "mutually_exclusive"

    return "alternative"


def _find_alternative_paths(graph,
                            sources,
                            destinies,
                            canonical_path,
                            min_genes=1):
    """
    Return a list of alternative paths.
    """
    paths = {'CanonicalPath': [], 'AlternativePath': [], 'ASE': []}
    canonical_set = set(canonical_path)
    for (src, i) in sources.items():
        for (dst, j) in destinies.items():
            if i < j:
                for path in nx.all_simple_paths(graph, src, dst):
                    if _is_alternative_path(graph, path, canonical_set, i, j,
                                            min_genes):
                        paths['CanonicalPath'].append(",".join(
                            canonical_path[i:j + 1]))
                        paths['AlternativePath'].append(",".join(path))
                        paths['ASE'].append(
                            _define_aes(graph, canonical_path, i, j, path))
    return pd.DataFrame(paths)


def conserved_ases(table, graph_file_name, min_genes=1):
    """
    Return a DataFrame of the conserved alternative splicing events detected.
    """
    graph = nx.read_gml(graph_file_name)
    canonical_path = get_canonical_path(table, graph)
    sources, destinies = _get_sources_and_destinies(graph, canonical_path,
                                                    min_genes)
    return _find_alternative_paths(graph, sources, destinies, canonical_path,
                                   min_genes)
