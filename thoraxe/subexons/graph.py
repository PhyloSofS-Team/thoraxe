"""
graph: Functions to store the splice graph with conservation information.
"""

import collections
import logging

from functools import reduce


def _gene2n_transcripts(data):
    """
    Return a dictionary from gene id to number of transcripts in the gene.
    """
    return data.groupby(
        'GeneID').apply(lambda df: len(set(df.TranscriptIDCluster))).to_dict()


def _edge2stats(edge2gene_list):
    """
    Return two dictionaries.

    - edge to set of genes
    - edge to dict from gene to number of transcripts in the gene with that edge
    """
    edge2genes = {}
    edge2gene_abundance = {}  # edge to transcript abundance in each gene

    for edge, gene_list in edge2gene_list.items():
        edge2genes[edge] = set(gene_list)
        edge2gene_abundance[edge] = collections.Counter(gene_list)

    return edge2genes, edge2gene_abundance


def _transcript_weighted_conservation(gene2transcriptnumber,
                                      edge2gene_abundance):
    """"
    Return the transcript weighted conservation score of each edge as a dict.
    """
    edge2trx_cons = {}
    n_genes = len(gene2transcriptnumber)
    for (edge, gene_abundance) in edge2gene_abundance.items():
        transcript_fraction_sum = 0.0
        for (gene, transcript_abundance) in gene_abundance.items():
            transcript_fraction_sum += (transcript_abundance /
                                        gene2transcriptnumber[gene])
        edge2trx_cons[edge] = transcript_fraction_sum / n_genes

    return edge2trx_cons


def nodes_and_edges2genes_and_transcripts(  # pylint: disable=too-many-locals
        data):
    """
    Return five dictionaries:

    - nodes to genes
    - edges to genes
    - nodes to transcripts
    - edges to transcripts
    - edges to the transcript weighted conservation score
    """
    gene2transcriptnumber = _gene2n_transcripts(data)

    node2genes = collections.defaultdict(set)
    edge2gene_list = collections.defaultdict(list)
    node2transcripts = collections.defaultdict(set)
    edge2transcripts = collections.defaultdict(set)
    for gene_id, gene in data.groupby('GeneID'):  # noqa pylint: disable=too-many-nested-blocks
        for transcript_id, transcript in gene.groupby('TranscriptIDCluster'):
            subexons = transcript.sort_values('SubexonRank')['S_exons']

            previous = 'start'
            node2genes[previous].update({gene_id})
            node2transcripts[previous].update({transcript_id})
            for _, subexon in enumerate(subexons):
                orthologs = subexon.split('/')
                for ortholog in orthologs:
                    node2genes[ortholog].update({gene_id})
                    edge2gene_list[(previous, ortholog)].append(gene_id)
                    node2transcripts[ortholog].update({transcript_id})
                    edge2transcripts[(previous,
                                      ortholog)].update({transcript_id})
                    previous = ortholog

            node2genes['stop'].update({gene_id})
            edge2gene_list[(previous, 'stop')].append(gene_id)
            node2transcripts['stop'].update({transcript_id})
            edge2transcripts[(previous, 'stop')].update({transcript_id})

    edge2genes, edge2gene_abundance = _edge2stats(edge2gene_list)

    edge2trx_cons = _transcript_weighted_conservation(gene2transcriptnumber,
                                                      edge2gene_abundance)

    return (node2genes, edge2genes, node2transcripts, edge2transcripts,
            edge2trx_cons)


def _get_elements(node2elements):
    """
    Return a set of elements from the node to (set of) elements dictionary.

    >>> sorted(_get_elements({'a':{1, 2}, 'b':{2, 3}, 'c':{1}, 'd':{4}}))
    [1, 2, 3, 4]
    """
    return reduce(set.union, node2elements.values())


def splice_graph_gml(  # pylint: disable=too-many-locals, too-many-arguments
        filename, node2genes, edge2genes, node2transcripts, edge2transcripts,
        edge2trx_cons, s_exon_2_char):
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

    n_genes = len(_get_elements(node2genes))
    n_transcripts = len(_get_elements(node2transcripts))

    with open(filename, 'w') as gml:
        gml.write('''
        graph [
            directed 1
            id 42
            label "splice graph of s-exons"
        ''')
        node2id = {}
        node_id = 1
        for node, genes in node2genes.items():
            conservation = 100.0 * (len(genes) / n_genes)
            transcripts = node2transcripts[node]
            transcript_fraction = 100.0 * (len(transcripts) / n_transcripts)
            genes_str = ','.join(sorted(genes))
            transcripts_str = ','.join(sorted(transcripts))
            gml.write('''
                node [
                    id {}
                    label "{}"
                    transcript_fraction {}
                    conservation {}
                    genes "{}"
                    transcripts "{}"
            '''.format(node_id, node, transcript_fraction, conservation,
                       genes_str, transcripts_str))
            if s_exon_2_char and node not in ['start', 'stop']:
                gml.write('''
                        phylosofs "{}"
                '''.format(s_exon_2_char[node]))
            gml.write('''
                ]
            ''')
            node2id[node] = node_id
            node_id += 1
        for edge, genes in edge2genes.items():
            conservation = 100.0 * (len(genes) / n_genes)
            transcripts = edge2transcripts[edge]
            transcript_fraction = 100.0 * (len(transcripts) / n_transcripts)
            transcript_weighted_conservation = edge2trx_cons[edge]
            genes_str = ','.join(sorted(genes))
            transcripts_str = ','.join(sorted(transcripts))
            gml.write('''
                edge [
                    source {}
                    target {}
                    transcript_fraction {}
                    conservation {}
                    transcript_weighted_conservation {}
                    genes "{}"
                    transcripts "{}"
                ]
            '''.format(node2id[edge[0]], node2id[edge[1]], transcript_fraction,
                       conservation, transcript_weighted_conservation,
                       genes_str, transcripts_str))
        gml.write('''
        ]
        ''')

    return filename
