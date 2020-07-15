"""
ases: Function to detect conserved alternative splicing events.

It finds the canonical path in the splice graph to detect conserved ASEs.
"""

import collections
import os

import pandas as pd
import networkx as nx
import numpy as np

from thoraxe.subexons import alignment

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
        s_exon_table,
        transcript_table,
        graph,
        column_order=None,
        delim='/'):
    """
    Return a DataFrame with the needed data to choose the canonical path.
    """
    data = {
        'GeneID': [],
        'TranscriptIDCluster': [],
        'TranscriptLength': [],
        'MinimumTranscriptWeightedConservation': [],
        'MinimumTranscriptFraction': [],
        'MinimumConservation': [],
        'TranscriptWeightedConservationSum': [],
        'TranscriptFractionSum': [],
        'ConservationSum': [],
        'MeanTranscriptWeightedConservation': [],
        'MeanTranscriptFraction': [],
        'MeanConservation': [],
        'TSL': [],
        'Path': []
    }

    transcript2tsl = transcript_table.loc[:,
                                          ["TranscriptIDCluster", "TSL"
                                           ]].set_index("TranscriptIDCluster"
                                                        ).to_dict()['TSL']

    for (trx, subdf) in s_exon_table.groupby('TranscriptIDCluster'):
        n_rows = len(subdf)
        s_exon_len = [len(_get_seq(subdf.S_exon_Sequence.iloc[0]))]
        path = ["start", subdf.S_exonID.iloc[0]]
        score = [
            graph.get_edge_data(
                'start',
                subdf.S_exonID.iloc[0])['transcript_weighted_conservation']
        ]
        conservation = [
            graph.get_edge_data('start',
                                subdf.S_exonID.iloc[0])['conservation']
        ]
        transcript_fraction = [
            graph.get_edge_data('start',
                                subdf.S_exonID.iloc[0])['transcript_fraction']
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
                conservation.append(
                    graph.get_edge_data(s_exon_1, s_exon_2)['conservation'])
                transcript_fraction.append(
                    graph.get_edge_data(s_exon_1,
                                        s_exon_2)['transcript_fraction'])

        score.append(
            graph.get_edge_data(subdf.S_exonID.iloc[n_rows - 1],
                                'stop')['transcript_weighted_conservation'])
        conservation.append(
            graph.get_edge_data(subdf.S_exonID.iloc[n_rows - 1],
                                'stop')['conservation'])
        transcript_fraction.append(
            graph.get_edge_data(subdf.S_exonID.iloc[n_rows - 1],
                                'stop')['transcript_fraction'])
        path.append('stop')

        data['GeneID'].append(subdf.GeneID.iloc[0])
        data['TranscriptIDCluster'].append(trx)
        data['TranscriptLength'].append(sum(s_exon_len))
        data['MinimumTranscriptWeightedConservation'].append(min(score))
        data['MinimumConservation'].append(min(conservation))
        data['MinimumTranscriptFraction'].append(min(transcript_fraction))
        data['TranscriptWeightedConservationSum'].append(sum(score))
        data['ConservationSum'].append(sum(conservation))
        data['TranscriptFractionSum'].append(sum(transcript_fraction))
        data['MeanTranscriptWeightedConservation'].append(np.mean(score))
        data['MeanConservation'].append(np.mean(conservation))
        data['MeanTranscriptFraction'].append(np.mean(transcript_fraction))
        data['TSL'].append(transcript2tsl[trx])
        data['Path'].append(delim.join(path))

    data_frame = pd.DataFrame(data)
    path2ngenes = data_frame.groupby('Path').apply(
        lambda df: len(set(df.GeneID))).to_dict()
    data_frame = data_frame.assign(
        PathGeneNumber=[path2ngenes[path] for path in data_frame.Path])

    if column_order is None:
        column_order = [
            'MinimumConservation', 'MinimumTranscriptWeightedConservation',
            'MeanTranscriptWeightedConservation', 'TranscriptLength', 'TSL'
        ]

    data_frame.drop_duplicates(inplace=True)
    data_frame.sort_values(column_order, ascending=False, inplace=True)

    return data_frame.reindex(columns=[
        'GeneID', 'TranscriptIDCluster', 'TranscriptLength', 'TSL', 'Path',
        'MinimumConservation', 'MinimumTranscriptFraction',
        'MinimumTranscriptWeightedConservation', 'ConservationSum',
        'TranscriptFractionSum', 'TranscriptWeightedConservationSum',
        'MeanConservation', 'MeanTranscriptFraction',
        'MeanTranscriptWeightedConservation', 'PathGeneNumber'
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


def _are_eclusives(  # pylint: disable=too-many-arguments
        s_exon_a, s_exon_b, path_list, path_dict, min_genes, min_transcripts,
        exclusives, s_exon2path):
    """
    Return True if both s-exons are mutually exclusive.

    `exclusives` and `s_exon2path` are dictionaries used to memoize results.
    """
    key = (s_exon_a, s_exon_b)
    if key not in exclusives:
        a_paths = _get_s_exon_paths(s_exon_a, path_list, s_exon2path)
        b_paths = _get_s_exon_paths(s_exon_b, path_list, s_exon2path)
        common_paths = a_paths.intersection(b_paths)
        paths_data = _genes_and_transcripts(common_paths, path_dict)
        exclusives[key] = not _min_filter(paths_data, min_genes,
                                          min_transcripts)

    return exclusives[key]


def _exclusive(  # pylint: disable=too-many-arguments
        canonical_s_exons, alternative_s_exons, path_list, path_dict,
        min_genes, min_transcripts, exclusives, s_exon2path):
    """
    Return a three elements tuple about the mutually exclusivity of the ASE.

    `exclusives` and `s_exon2path` are dictionaries used to memoize results.
    """
    exclusive_canonical = []
    exclusive_alternative = []
    me_status = "mutually_exclusive"
    for alt in alternative_s_exons:
        for can in canonical_s_exons:
            if _are_eclusives(can, alt, path_list, path_dict, min_genes,
                              min_transcripts, exclusives, s_exon2path):
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


def _genes_and_transcripts(paths, path_dict):
    """
    Return a PathData object with the genes and transcripts in `paths`.
    """
    genes = set()
    transcripts = set()
    for path in paths:
        path_data = path_dict[path]
        genes.update(path_data.Genes)
        transcripts.update(path_data.Transcripts)
    return PathData(genes, transcripts)


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
    n_genes = len(gene2transcripts)
    transcript_fraction_sum = 0.0
    for gene in subpath_data.Genes:
        gene_transcripts = gene2transcripts[gene]
        subpath_transcripts = len(
            gene_transcripts.intersection(subpath_data.Transcripts))
        gene_transcripts = len(gene_transcripts)
        transcript_fraction_sum += (subpath_transcripts / gene_transcripts)
    return transcript_fraction_sum / n_genes


def _min_filter(path_data, min_genes, min_transcripts):
    """
    Return `True` is path_data has more than min_genes and min_transcripts.
    """
    if min_genes > 1:
        if len(path_data.Genes) < min_genes:
            return False
    if min_transcripts > 1:
        if len(path_data.Transcripts) < min_transcripts:
            return False
    return True


def detect_ases(  # pylint: disable=too-many-locals,too-many-nested-blocks
        path_table,
        min_genes=0,
        min_transcripts=1,
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
                if _min_filter(alternative_data, min_genes, min_transcripts):
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
                                path_list, path_dict, min_genes,
                                min_transcripts, exclusives, s_exon2path)
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


def _get_path_sequence(s_exon_df, path, path_genes, delim='/'):
    """
    Return the residue sequences for the path.
    """
    s_exons = path.split(delim)
    sequences = []
    for gene in path_genes:
        gene_rows = s_exon_df.GeneID == gene
        sequence = []
        for s_exon in s_exons:
            seqs = s_exon_df.loc[(s_exon_df.S_exonID == s_exon) & gene_rows,
                                 'S_exon_Sequence'].drop_duplicates()
            if seqs.size:
                sequence.append(seqs.iloc[0])
        sequences.append(''.join(sequence))
    return sequences


def _len_stats(sequences):
    """
    Return the minimum, mean and maximum sequence length.

    >>> _len_stats(['A', 'AA', 'AAA'])
    (1, 2.0, 3)
    """
    lengths = [len(seq) for seq in sequences]
    return min(lengths), sum(lengths) / len(lengths), max(lengths)


def _get_path_consensus(s_exon_msas, path, path_genes, delim='/'):
    """
    Return the consensus sequence for the given path in the given genes.
    """
    consensus = []
    for s_exon in path.split(delim):
        if s_exon in s_exon_msas:
            msa = s_exon_msas[s_exon]
            sub_msa = alignment.get_submsa(msa, path_genes)
            consensus.append(alignment.get_consensus(sub_msa))
        else:
            consensus.append('')
    return ''.join(consensus)


def create_ases_table(s_exon_msas, s_exon_table, events, delim='/'):
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
        'CanonicalPathMinSeqLength': [],
        'CanonicalPathMeanSeqLength': [],
        'CanonicalPathMaxSeqLength': [],
        'AlternativePathMinSeqLength': [],
        'AlternativePathMeanSeqLength': [],
        'AlternativePathMaxSeqLength': [],
        'CanonicalPathConsensus': [],
        'AlternativePathConsensus': [],
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
        can_stats = _len_stats(
            _get_path_sequence(s_exon_table, key[0], can_genes, delim=delim))
        paths['CanonicalPathMinSeqLength'].append(can_stats[0])
        paths['CanonicalPathMeanSeqLength'].append(can_stats[1])
        paths['CanonicalPathMaxSeqLength'].append(can_stats[2])
        alt_stats = _len_stats(
            _get_path_sequence(s_exon_table, key[1], alt_genes, delim=delim))
        paths['AlternativePathMinSeqLength'].append(alt_stats[0])
        paths['AlternativePathMeanSeqLength'].append(alt_stats[1])
        paths['AlternativePathMaxSeqLength'].append(alt_stats[2])
        paths['CanonicalPathConsensus'].append(
            _get_path_consensus(s_exon_msas, key[0], can_genes, delim=delim))
        paths['AlternativePathConsensus'].append(
            _get_path_consensus(s_exon_msas, key[1], alt_genes, delim=delim))
    data_frame = pd.DataFrame(paths)
    data_frame.sort_values(
        ['NumberOfCommonGenes', 'CanonicalPath', 'AlternativePath'],
        ascending=False,
        inplace=True)
    data_frame.drop(columns=['NumberOfCommonGenes'], inplace=True)
    return data_frame


def conserved_ases(  # pylint: disable=too-many-arguments
        s_exon_table,
        transcript_table,
        s_exon_msas,
        graph_file_name,
        min_genes=1,
        min_transcripts=2,
        column_order=None,
        delim='/'):
    """
    Return two DataFrames.

    - A table of transcripts/paths.
    - A table of the conserved alternative splicing events detected.
    """
    graph = read_splice_graph(graph_file_name)
    path_table = get_transcript_scores(s_exon_table,
                                       transcript_table,
                                       graph,
                                       column_order=column_order)
    events = detect_ases(path_table,
                         min_genes=min_genes,
                         min_transcripts=min_transcripts,
                         delim=delim)
    ases_table = create_ases_table(s_exon_msas,
                                   s_exon_table,
                                   events,
                                   delim=delim)

    return path_table, ases_table
