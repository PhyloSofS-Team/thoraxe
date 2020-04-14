from collections import namedtuple
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
        ase_type = "insertion"
    elif len(ase_alternative) == 2:
        ase_type = "deletion"
    elif ase_alternative[0] == "start":
        if ase_alternative[-1] == "stop":
            ase_type = "fully_alternative"
        else:
            ase_type = "alternative_start"
    elif ase_alternative[-1] == "stop":
        ase_type = "alternative_end"
    else:
        ase_type = "alternative"
    return(ase_type)


def _get_s_exon_paths(s_exon, path_list, s_exon2path):
    """
    Return the set of path with a given s_exon.

    It uses s_exon2path to memoize the result.
    """
    if s_exon not in s_exon2path:
        s_exon2path[s_exon] = set(filter(lambda x: '/' + s_exon + '/' in x, path_list)))

    return s_exon2path[s_exon]


def _are_eclusives(s_exon_a, s_exon_b, path_list, exclusives, s_exon2path):
    key = (s_exon_a, s_exon_b)
    if key not in exclusives:
        a_paths = _get_s_exon_paths(s_exon_a, path_list, s_exon2path)
        b_paths = _get_s_exon_paths(s_exon_b, path_list, s_exon2path)
        common_paths = a_paths.intersection(b_paths)
        exclusives[key] = len(common_paths) == 0

    return exclusives[key]


def _exclusive(canonical_s_exons, alternative_s_exons, path_list, exclusives, s_exon2path):
    exclusive_canonical = []
    exclusive_alternative = []
    me_status = "mutually_exclusive"
    for alt in alternative_s_exons:
        for can in canonical_s_exons:
            if _are_eclusives(s_exon_a, s_exon_b, path_list, exclusives, s_exon2path):
                exclusive_canonical.append(can)
                exclusive_alternative.append(alt)
            else:
                me_status = "partially_mutually_exclusive"

    return exclusive_canonical, exclusive_alternative, me_status


def _exclusive(paths, canonical_s_exons, alternative_s_exons):
    selected_paths = set()
    # we initialize the status at partially mutually exclusive
    me = "pMe"
    for se in canonical_s_exons:
        selected_paths = selected_paths.union(
            set(filter(lambda x: '/' + se + '/' in x, paths)))
    alternative_filtered = []
    for se in alternative_s_exons:
        if not list(filter(lambda x: '/' + se + '/' in x, selected_paths)):
            alternative_filtered.append(se)
    if len(alternative_filtered) > 0:
        # if none of the elements in se were removed
        if len(alternative_s_exons) == len(alternative_filtered):
            # then the sets are fully mutually exclusive
            me = "ME"
        # conditional return: if the alternative_filtered is empty, the function returns none
        return(','.join(canonical_s_exons), ','.join(alternative_filtered), me)


PathData = namedtuple('PathData', ['Genes', 'Transcripts'])


def _get_path_dict(path_table):
    """
    Return a list of paths and a Dict from a path to its PathData namedtuple.

    Paths are reperesented with strings, where each node is separated by `/`.
    The path list is useful to keep path order and to make paths accessible by
    index. Each `PathData` `namedtuple` in the `dict` contains the list of
    genes and `TranscriptIDCluster`s.
    """
    path_list = []
    path_dict = dict()
    for row in path_table.itertuples():
        if row.Path not in path_dict:
            path_list.append(row.Path)
            path_dict[row.Path] = PathData(
                {row.GeneID}, {row.TranscriptIDCluster})
        else:
            path_dict[row.Path].Genes.add(row.GeneID)
            path_dict[row.Path].Transcripts.add(row.TranscriptIDCluster)
    return path_list, path_dict


def _subpath_data(subpath, path_dict, subpath_dict):
    """
    Return the `PathData` for the given subpath.

    This function memoizes using `subpath_dict`.
    """
       if subpath not in subpath_dict:
            subpath_dict[subpath] = PathData(set(), set())
            for (path, path_data) in path_dict.items():
                if subpath in path:
                    subpath_dict[subpath].Genes.update(path_data.Genes)
                    subpath_dict[subpath].Transcripts.update(
                        path_data.Transcripts)
        return subpath_dict[subpath]


def _detect_ase(path_table, min_genes=1, delim='/'):
    subpath_dict = {}  # memoization dict for _subpath_data
    exclusives = {}  # memoization dict for _exclusive
    s_exon2path = {}  # memoization dict for _exclusive
    path_list, path_dict = _get_path_dict(path_table)
    # the first path is the canonical path
    canonical = path_list[0].split('/')
    canonical_node2index = {node : idx for (idx, node) in enumerate(canonical)}
    # create a dico of events, where each key is a tuple
    # (canonical_path, alternative_path) and the value is a tuple of:
    # - a PathData tuple for each path
    # - the type of alternative splicing event
    # - a note about the mutually exclusivity of the event
    # - the mutually exclusive s-exons in the canonical
    # - and in the alternative paths
    events = {}
    for joined_path in path_list[1:]:
        possibly_exclusive = set(["alternative_start", "alternative", "alternative_end", "fully_alternative"])
        path = joined_path.split("/")
        # we start at 1 because 0 is necessarily the "start"
        i = 1
        j = 1
        path_len = len(path)
        while i < path_len:
            if path[i] != canonical[j]:
                istop = i
                jstop = j
                # find the next common s-exons
                for k in range(istop, path_len):
                    can_index = canonical_node2index.get(path[k], -1)
                    istop = k
                    if can_index >= jstop:
                        jstop = can_index  # its index in the canonical path
                        break
                # create alternative path
                ase_alternative = path[(i - 1):(istop + 1)]
                joined_alternative = "/".join(ase_alternative)
                # check whether all its s-exons are present in at least one
                # species
                alternative_data = _subpath_data(joined_alternative, path_dict,
                    subpath_dict)
                if len(alternative_data.Genes) > min_genes:
                    # create the canonical (main) path
                    ase_canonical = canonical[(j - 1):(jstop + 1)]
                    joined_canonical = "/".join(ase_canonical)
                    # create the event as a tuple of canonical (main) and
                    # alternative paths
                    event = (joined_canonical, joined_alternative)
                    # create it in the dictionary if it is detected for the
                    # first time
                    if event not in events:
                        # get data for the canonical path
                        canonical_data = _subpath_data(joined_canonical,
                            path_dict, subpath_dict)
                        # determine the ASE type according the classical
                        # classification
                        ase_type = _ase_type(ase_canonical, ase_alternative)
                        # determine whether part of all of the alternative path
                        # is mutually exclusive with the canonical (main) path
                        me_status = ""
                        if ase_type in possibly_exclusive:
                            canonical_s_exons = ase_canonical[1:-1]
                            alternative_s_exons = ase_alternative[1:-1]
                            exclusive_can, exclusive_alt, me_status = _exclusive(
                                canonical_s_exons, alternative_s_exons,
                                path_list, exclusives, s_exon2path)
                        # TO DO
                        events[event] = (alternative_data, canonical_data,
                            me_status, exclusive_can, exclusive_alt)
                i = istop + 1
                j = jstop + 1
            else:
                i += 1
                j += 1
    return events


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
