"""
ases: Function to detect conserved alternative splicing events.

It finds the canonical path in the splice graph to detect conserved ASEs.
"""

import collections
import networkx as nx
import pandas as pd

PathData = collections.namedtuple('PathData', ['Genes', 'Transcripts'])


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
    path2ngenes = data_frame.groupby('Path').apply(
        lambda df: len(set(df.GeneID))).to_dict()
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
    """
    Return the type of ASE, it does not define mutual exclusivity.
    """
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
    return ase_type


def _get_s_exon_paths(s_exon, path_list, s_exon2path):
    """
    Return the set of path with a given s_exon.

    It uses s_exon2path to memoize the result.
    """
    if s_exon not in s_exon2path:
        s_exon2path[s_exon] = set(
            filter(lambda x: '/' + s_exon + '/' in x, path_list))

    return s_exon2path[s_exon]


def _are_eclusives(s_exon_a, s_exon_b, path_list, exclusives, s_exon2path):
    """
    Return True if both s-exons are mutually exclusive.

    `exclusives` and `s_exon2path` are dictionaries used to memoize results.
    """
    key = (s_exon_a, s_exon_b)
    if key not in exclusives:
        a_paths = _get_s_exon_paths(s_exon_a, path_list, s_exon2path)
        b_paths = _get_s_exon_paths(s_exon_b, path_list, s_exon2path)
        common_paths = a_paths.intersection(b_paths)
        exclusives[key] = len(common_paths) == 0

    return exclusives[key]


def _exclusive(canonical_s_exons, alternative_s_exons, path_list, exclusives,
               s_exon2path):
    """
    Return a three elements tuple about the mutually exclusivity of the ASE.
    """
    exclusive_canonical = []
    exclusive_alternative = []
    me_status = "mutually_exclusive"
    for alt in alternative_s_exons:
        for can in canonical_s_exons:
            if _are_eclusives(can, alt, path_list, exclusives, s_exon2path):
                exclusive_canonical.append(can)
                exclusive_alternative.append(alt)
            else:
                me_status = "partially_mutually_exclusive"
    if len(exclusive_canonical) == 0:  # == exclusive_alternative
        me_status = ""

    return exclusive_canonical, exclusive_alternative, me_status


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
            path_dict[row.Path] = PathData({row.GeneID},
                                           {row.TranscriptIDCluster})
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
                subpath_dict[subpath].Transcripts.update(path_data.Transcripts)
    return subpath_dict[subpath]


def _gene2transcripts(path_table):
    """
    Return a Dict from `GeneID` to `TranscriptIDCluster` from the `path_table`.
    """
    gene2transcripts = collections.defaultdict(set)
    for row in path_table.itertuples():
        gene2transcripts[row.GeneID].add(row.TranscriptIDCluster)
    return gene2transcripts


def _transcript_weighted_conservation(subpath_data, gene2transcripts):
    """
    Return the transcript weighted conservation of a given subpath.
    """
    acc = 0.0
    for gene in subpath_data.Genes:
        gene_transcripts = gene2transcripts[gene]
        acc += (
            (len(gene_transcripts.intersection(subpath_data.Transcripts))) /
            len(gene_transcripts))
    return acc / len(gene2transcripts)


def detect_ases(  # pylint: disable=too-many-locals,too-many-nested-blocks
        path_table,
        min_genes=1,
        delim='/'):
    """
    Return a dictionary from ASEs to their data.

    It takes as imput a pandas DataFrame object with the path table.
    Each key of the output dictionary is a tuple
    (canonical_path, alternative_path) and the value is a tuple of:

     - one PathData tuple for the canonical path
     - and one for the alternative path;
     - the transcript weighted conservation for the canonical sub-path,
     - and for the alternative sub-path,
     - the type of alternative splicing event, e.g. "insertion";
     - the mutually exclusivity of the event ("" for non mutually exclusive);
     - one list of mutually exclusive s-exons in canonical paths
     - and one for the alternative paths.

    The last two list are paired, i.e. the two first s-exons are mutually
    exclusive.
    """
    subpath_dict = {}  # memoization dict for _subpath_data
    exclusives = {}  # memoization dict for _exclusive
    s_exon2path = {}  # memoization dict for _exclusive
    gene2transcripts = _gene2transcripts(path_table)
    path_list, path_dict = _get_path_dict(path_table)
    # the first path is the canonical path
    canonical = path_list[0].split(delim)
    canonical_node2index = {node: idx for (idx, node) in enumerate(canonical)}
    # create a dico of events, where each key is a tuple
    # (canonical_path, alternative_path) and the value is a tuple
    events = {}
    for joined_path in path_list[1:]:
        possibly_exclusive = set([
            "alternative_start", "alternative", "alternative_end",
            "fully_alternative"
        ])
        path = joined_path.split(delim)
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
                joined_alternative = delim.join(ase_alternative)
                # check whether all its s-exons are present in at least one
                # species
                alternative_data = _subpath_data(joined_alternative, path_dict,
                                                 subpath_dict)
                if len(alternative_data.Genes) > min_genes:
                    # create the canonical (main) path
                    ase_canonical = canonical[(j - 1):(jstop + 1)]
                    joined_canonical = delim.join(ase_canonical)
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
                        me_can = []
                        me_alt = []
                        if ase_type in possibly_exclusive:
                            canonical_s_exons = ase_canonical[1:-1]
                            alternative_s_exons = ase_alternative[1:-1]
                            me_can, me_alt, me_status = _exclusive(
                                canonical_s_exons, alternative_s_exons,
                                path_list, exclusives, s_exon2path)
                        # transcript weighted conservation for both sub-paths
                        twc_can = _transcript_weighted_conservation(
                            canonical_data, gene2transcripts)
                        twc_alt = _transcript_weighted_conservation(
                            alternative_data, gene2transcripts)
                        # store alternative splicing event
                        events[event] = (canonical_data, alternative_data,
                                         twc_can, twc_alt, ase_type, me_status,
                                         me_can, me_alt)
                i = istop + 1
                j = jstop + 1
            else:
                i += 1
                j += 1
    return events


def create_ases_table(events, delim='/'):
    """
    Create an ASE data frame from the dictionary returned by `detect_ases`.
    """
    paths = {
        'CanonicalPath': [],
        'AlternativePath': [],
        'ASE': [],
        'MutualExclusivity': [],
        'MutualExclusiveCanonical': [],
        'MutualExclusiveAlternative': [],
        'CanonicalPathTranscriptWeightedConservation': [],
        'AlternativePathTranscriptWeightedConservation': [],
        'CanonicalPathGenes': [],
        'AlternativePathGenes': [],
        'CommonGenes': [],
        'NumberOfCommonGenes': [],
    }
    for (key, value) in events.items():
        paths['CanonicalPath'].append(key[0])
        paths['AlternativePath'].append(key[1])
        paths['ASE'].append(value[4])
        paths['MutualExclusivity'].append(value[5])
        paths['MutualExclusiveCanonical'].append(delim.join(value[6]))
        paths['MutualExclusiveAlternative'].append(delim.join(value[7]))
        paths['CanonicalPathTranscriptWeightedConservation'].append(value[2])
        paths['AlternativePathTranscriptWeightedConservation'].append(value[3])
        can_genes = value[0].Genes
        alt_genes = value[1].Genes
        paths['CanonicalPathGenes'].append(delim.join(sorted(can_genes)))
        paths['AlternativePathGenes'].append(delim.join(sorted(alt_genes)))
        common_genes = can_genes.intersection(alt_genes)
        paths['CommonGenes'].append(delim.join(sorted(common_genes)))
        paths['NumberOfCommonGenes'].append(len(common_genes))
    data_frame = pd.DataFrame(paths)
    data_frame.sort_values(
        ['NumberOfCommonGenes', 'CanonicalPath', 'AlternativePath'],
        ascending=False,
        inplace=True)
    data_frame.drop(columns=['NumberOfCommonGenes'], inplace=True)
    return data_frame


def conserved_ases(table, graph_file_name, min_genes=1, delim='/'):
    """
    Return two DataFrames.

    - A table of transcripts/paths.
    - A table of the conserved alternative splicing events detected.
    """
    graph = read_splice_graph(graph_file_name)
    path_table = get_transcript_scores(table, graph)
    events = detect_ases(path_table, min_genes=min_genes, delim=delim)
    ases_table = create_ases_table(events, delim=delim)
    return path_table, ases_table
