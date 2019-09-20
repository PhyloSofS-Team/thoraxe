"""
graph: Functions to store the splice graph with conservation information.
"""

import collections
import logging

from functools import reduce


def nodes_and_edges2genes(data):
    """
    Return two dictionaries, ones from nodes and other from edges to genes.
    """
    node2genes = collections.defaultdict(set)
    edge2genes = collections.defaultdict(set)
    for gene_id, gene in data.groupby('GeneID'):  # noqa pylint: disable=too-many-nested-blocks
        for _, transcript in gene.groupby('TranscriptID'):
            subexons = transcript.sort_values('SubexonRank')['HomologousExons']

            previous = 'start'
            node2genes[previous].update({gene_id})
            for _, subexon in enumerate(subexons):
                orthologs = subexon.split('/')
                for ortholog in orthologs:
                    node2genes[ortholog].update({gene_id})
                    edge2genes[(previous, ortholog)].update({gene_id})
                    previous = ortholog

            node2genes['stop'].update({gene_id})
            edge2genes[(previous, 'stop')].update({gene_id})

    return node2genes, edge2genes


def _get_genes(node2genes):
    """
    Return a set of genes from the node to (set of) genes dictionary.

    >>> sorted(_get_genes({'a':{1, 2}, 'b':{2, 3}, 'c':{1}, 'd':{4}}))
    [1, 2, 3, 4]
    """
    return reduce(set.union, node2genes.values())


def splice_graph_gml(filename, node2genes, edge2genes):
    """
    Store the splice graph in GML (Graph Modelling Language) format.

    It stores the conservation level information for nodes and edges, where
    conservation is the fraction of species that has that particular feature.
    """
    if not filename.endswith('.gml'):
        filename += '.gml'
        logging.warning(
            '.gml extension was added, the splice graph will be stored at %s',
            filename)

    n_genes = len(_get_genes(node2genes))

    with open(filename, 'w') as gml:
        gml.write('''
        graph [
            directed 1
            id 42
            label "splice graph of orthologous exon groups"
        ''')
        node2id = {}
        node_id = 1
        for node, genes in node2genes.items():
            conservation = 100.0 * (len(genes) / n_genes)
            gml.write('''
                node [
                    id {}
                    label "{}"
                    conservation {}
                ]
            '''.format(node_id, node, conservation))
            node2id[node] = node_id
            node_id += 1
        for edge, genes in edge2genes.items():
            conservation = 100.0 * (len(genes) / n_genes)
            gml.write('''
                edge [
                    source {}
                    target {}
                    conservation {}
                ]
            '''.format(node2id[edge[0]], node2id[edge[1]], conservation))
        gml.write('''
        ]
        ''')

    return filename
