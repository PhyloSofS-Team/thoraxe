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


def get_transcript_scores(  # pylint: disable=too-many-locals
        table, graph, delim='/'):
    """
    Return a DataFrame with the needed data to choose the canonical path.
    """
    bless_humans = 'homo_sapiens' in set(table.Species)
    data = {
        'GeneID': [],
        'TranscriptIDCluster': [],
        'TranscriptLength': [],
        'MinimumTranscriptWeightedConservation': [],
        'Path': []
    }
    if bless_humans:
        data['IsHuman'] = []

    for (trx, subdf) in table.groupby('TranscriptIDCluster'):
        n_rows = len(subdf)
        s_exon_len = [len(_get_seq(subdf.S_exon_Sequence.iloc[0]))]
        path = ["start", subdf.S_exonID.iloc[0]]
        score = [
            graph.get_edge_data(
                'start',
                subdf.S_exonID.iloc[0])['transcript_weighted_conservation']
        ]
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

        score.append(
            graph.get_edge_data(subdf.S_exonID.iloc[n_rows - 1],
                                'stop')['transcript_weighted_conservation'])
        path.append('stop')

        data['GeneID'].append(subdf.GeneID.iloc[0])
        data['TranscriptIDCluster'].append(trx)
        data['TranscriptLength'].append(sum(s_exon_len))
        data['MinimumTranscriptWeightedConservation'].append(min(score))
        data['Path'].append(delim.join(path))
        if bless_humans:
            data['IsHuman'].append(
                int(subdf.Species.iloc[0] == 'homo_sapiens'))

    data_frame = pd.DataFrame(data)
    path2ngenes = data_frame.groupby(
        'Path').apply(lambda df: len(set(df.GeneID))).to_dict()
    data_frame = data_frame.assign(
        PathGeneNumber=[path2ngenes[path] for path in data_frame.Path])

    column_order = [
        'PathGeneNumber', 'MinimumTranscriptWeightedConservation',
        'TranscriptLength'
    ]
    if bless_humans:
        column_order.insert(2, 'IsHuman')

    data_frame.drop_duplicates(inplace=True)
    data_frame.sort_values(column_order, ascending=False, inplace=True)

    if bless_humans:
        data_frame.drop(columns=['IsHuman'], inplace=True)

    return data_frame.reindex(columns=[
        'GeneID', 'TranscriptIDCluster', 'TranscriptLength', 'Path',
        'MinimumTranscriptWeightedConservation', 'PathGeneNumber'
    ])


def get_canonical_path(table, graph, delim='/'):
    """
    Return the path table and the canonical path.
    """
    trx_data = get_transcript_scores(table, graph)
    path = trx_data.Path.iloc[0]
    return trx_data, path.split(delim)


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


def _get_genes(graph, src, dst):
    """
    Return the set of genes in the edge (src, dst).
    """
    return set(graph.get_edge_data(src, dst)['genes'].split(','))


def _is_alternative_path(  # noqa pylint: disable=too-many-arguments
        graph,
        path,
        canonical_set,
        i,
        j,
        min_genes=1):
    """
    Return the set of genes if the given path is alternative to the canonical.

    Otherwise, it returns an empty set (False).
    """
    internal_nodes = path[1:-1]

    if (not internal_nodes) and (j - i == 1):
        return set()

    for node in internal_nodes:
        if node in canonical_set:
            return set()

    n_nodes = len(path)

    genes = _get_genes(graph, path[0], path[1])

    if len(genes) < min_genes:
        return set()

    for k in range(2, n_nodes):
        genes.intersection_update(_get_genes(graph, path[k - 1], path[k]))

    if len(genes) >= min_genes:
        return genes

    return set()


def _genes_canonical(memo, graph, canonical_path, i, j):
    """
    Return the set of genes in the canonical path, it does memoization.
    """
    if (i, j) not in memo:
        genes = _get_genes(graph, canonical_path[i], canonical_path[i + 1])
        for k in range(i + 1, j):
            genes.intersection_update(
                _get_genes(graph, canonical_path[k], canonical_path[k + 1]))
        memo[(i, j)] = genes

    return memo[(i, j)]


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


def _setup_graph(graph):
    """
    Return a copy of the graph with a distance atribute for all_shortest_paths.
    """
    new_graph = graph.copy(as_view=False)
    nx.set_edge_attributes(
        new_graph, {(s, t): {
            'distance':
            1.0
            - new_graph.get_edge_data(s, t)['transcript_weighted_conservation']
        }
            for (s, t) in new_graph.edges()})
    return new_graph


def _find_alternative_paths(  # noqa pylint: disable=too-many-arguments,too-many-locals
        splice_graph,
        sources,
        destinies,
        canonical_path,
        min_genes=1,
        delim='/'):
    """
    Return a list of alternative paths.
    """
    memo = {}
    paths = {
        'CanonicalPath': [],
        'AlternativePath': [],
        'ASE': [],
        'AlternativePathGeneNumber': [],
        'AlternativePathGenes': [],
        'CanonicalPathGeneNumber': [],
        'CanonicalPathGenes': []
    }
    canonical_set = set(canonical_path)
    graph = _setup_graph(splice_graph)
    for (src, i) in sources.items():
        for (dst, j) in destinies.items():
            if i < j:
                for path in nx.all_shortest_paths(graph,
                                                  src,
                                                  dst,
                                                  weight='distance'):
                    genes = _is_alternative_path(graph, path, canonical_set, i,
                                                 j, min_genes)
                    if genes:
                        canonical_genes = _genes_canonical(
                            memo, graph, canonical_path, i, j)
                        paths['CanonicalPath'].append(
                            delim.join(canonical_path[i:j + 1]))
                        paths['AlternativePath'].append(delim.join(path))
                        paths['ASE'].append(
                            _define_aes(graph, canonical_path, i, j, path))
                        paths['AlternativePathGeneNumber'].append(len(genes))
                        paths['AlternativePathGenes'].append(
                            delim.join(sorted(genes)))
                        paths['CanonicalPathGeneNumber'].append(
                            len(canonical_genes))
                        paths['CanonicalPathGenes'].append(
                            delim.join(sorted(canonical_genes)))
    data_frame = pd.DataFrame(paths)
    data_frame.sort_values(
        ['AlternativePathGeneNumber', 'CanonicalPath', 'AlternativePath'],
        ascending=False,
        inplace=True)
    return data_frame.reindex(columns=[
        'CanonicalPath', 'AlternativePath', 'ASE', 'CanonicalPathGeneNumber',
        'CanonicalPathGenes', 'AlternativePathGeneNumber',
        'AlternativePathGenes'
    ])


def read_splice_graph(graph_file_name):
    """
    Read ThorAxe's splice graph using NetworkX.
    """
    with open(graph_file_name, 'r', encoding="utf-8") as gml:
        lines = [
            line.encode('ascii') for line in gml if "phylosofs" not in line
        ]
    return nx.read_gml(lines)


def _ase_type(ase_canonical, ase_alternative):
    ase_type = ""
    if len(ase_canonical) == 2:
        ase_type = "Insert"
    elif len(ase_alternative) == 2:
        ase_type = "Del"
    elif ase_alternative[0] == "start":
        if ase_alternative[-1] == "stop":
            ase_type = "Full"
        else:
            ase_type = "AlterS"
    elif ase_alternative[-1] == "stop":
        ase_type = "AlterE"
    else:
        ase_type = "AlterI"
    return(ase_type)


def _exclusive(paths, set1, set2):
    selected_paths = set()
    # we initialize the status at partially muually exclusive
    me = "pMe"
    for se in set1:
        selected_paths = selected_paths.union(
            set(filter(lambda x: '/' + se + '/' in x, paths)))
    set2filt = []
    for se in set2:
        if not list(filter(lambda x: '/' + se + '/' in x, selected_paths)):
            set2filt.append(se)
    if len(set2filt) > 0:
        # if none of the elements in se were removed
        if len(set2) == len(set2filt):
            # then the sets are fully mutually exclusive
            me = "ME"
        # conditional return: if the set2filt is empty, the function returns none
        return(','.join(set1), ','.join(set2filt), me)


def _detect_ase(path_table):
    paths = path_table.Path
    canonical = paths[0].split('/')
    # create a dico of events, where each key is a tuple (canonical_path, alternative_path)
    # and the values are the number of occurrences
    events = {}
    exclu = set()
    for joined_path in paths[1:]:
        path = joined_path.split("/")
        # we start at 1 because 0 is necessarily the start
        i = 1
        j = 1
        while i < len(path):
            if path[i] != canonical[j]:
                istop = i
                jstop = j
                # find the next common s-exons
                while not path[istop] in canonical[jstop:]:
                    istop += 1
                # find its index in the canonical path
                jstop = canonical.index(path[istop])
                # create alternative path
                ase_alternative = path[(i - 1):(istop + 1)]
                # check whether all its s-exons are present in at least one species
                # check s-exon freq here! # set(ase_alternative).issubset(de.keys()): # set to True if you want all events
                if True:
                    # create the canonical (main) path
                    ase_canonical = canonical[(j - 1):(jstop + 1)]
                    # create the event as a tuple of canonical (main) and alternative paths
                    event = ("/".join(ase_canonical),
                             "/".join(ase_alternative))
                    # max number of s-exons involved
                    ne = max(len(ase_canonical[1:-1]),
                             len(ase_alternative[1:-1]))
                    # update the number of occurrences of the event
                    if event in events:
                        events[event][0] += 1
                    # create it in the dictionary if it is detected for the first time
                    else:
                        # determine its type according the classical classification
                        ase_type = _ase_type(ase_canonical, ase_alternative)
                        # number of occurrences, type of ase, mutually exclusive status
                        # the latter will remain "-" in case of deletion, insertion, and alternative but not exclusive at all
                        events[event] = [1, ase_type, "-", ne, 0]
                        # determine whether part of all of the alternative path is
                        # mutually exclusive with the canonical (main) path
                        if ase_type in ["AlterS", "AlterI", "AlterE"]:
                            # ratio between the max and the other one
                            events[event][4] = max(
                                ne / len(ase_canonical[1:-1]), ne / len(ase_alternative[1:-1]))
                            aseExclu = _exclusive(
                                paths, ase_canonical[1:-1], ase_alternative[1:-1])
                            if aseExclu:
                                exclu.add(aseExclu[:2])
                                events[event][2] = aseExclu[2]
                i = istop + 1
                j = jstop + 1
            else:
                i += 1
                j += 1
    return (events, exclu)


def conserved_ases(table, graph_file_name, min_genes=1, delim='/'):
    """
    Return two DataFrames.

    - A table of transcripts/paths.
    - A table of the conserved alternative splicing events detected.
    """
    graph = read_splice_graph(graph_file_name)
    path_table, canonical_path = get_canonical_path(table, graph)
    sources, destinies = _get_sources_and_destinies(graph, canonical_path,
                                                    min_genes)
    return path_table, _find_alternative_paths(graph, sources, destinies,
                                               canonical_path, min_genes,
                                               delim)
