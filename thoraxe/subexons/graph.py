"""
graph: Functions to store the splice graph with conservation information.
"""

import collections


def nodes_and_edges2genes(data):
    """
    Return two dictionaries, ones from nodes and other from edges to genes.
    """
    node2genes = collections.defaultdict(set)
    edge2genes = collections.defaultdict(set)
    for gene_id, gene in data.groupby('Gene stable ID'):  # noqa pylint: disable=too-many-nested-blocks
        for _, transcript in gene.groupby('Transcript stable ID'):
            subexons = transcript.sort_values(
                'Subexon rank in transcript')['HomologousExons']
            previous = ''
            for _, subexon in enumerate(subexons):
                orthologs = subexon.split('-')
                for j, ortholog in enumerate(orthologs):
                    if j == 0:
                        if previous != '':
                            node2genes[previous].update({gene_id})
                            node2genes[ortholog].update({gene_id})
                            edge2genes[(previous, ortholog)].update({gene_id})
                        previous = ortholog
                        continue
                    node2genes[previous].update({gene_id})
                    node2genes[ortholog].update({gene_id})
                    edge2genes[(previous, ortholog)].update({gene_id})
                    previous = ortholog
    return node2genes, edge2genes
