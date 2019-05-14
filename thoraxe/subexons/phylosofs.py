"""
phylosofs: Functions to generate PhyloSofS input.
"""

import string

import pandas as pd

from Bio import Phylo

CHARS = [
    char for char in string.printable
    if char not in {' ', '\t', '\n', '\r', '\x0b', '\x0c', '\\'}
]

MAX_EXONS = len(CHARS)


def get_exon2char(exons):
    """
    Return a dictionary from orthologous exon group to a single character.

    >>> result = get_exon2char(['1_0-1_1', '2_0'])
    >>> sorted(result)
    ['1_0', '1_1', '2_0']
    >>> [result[exon] for exon in sorted(result)]
    ['0', '1', '2']
    """
    homologous_exons = {}
    i = 0
    for exon in exons:
        for subexon in str(exon).split('-'):
            if subexon not in homologous_exons:
                if i < MAX_EXONS:
                    homologous_exons[subexon] = CHARS[i]
                    i += 1
                else:
                    raise Exception(
                        'PhyloSofS can parse more than {} homologous exons.'.
                        format(MAX_EXONS))
    return homologous_exons


def get_transcript2phylosofs(data, exon2char):
    """
    Return a dictionary from transcript id to its phylosofs representation.
    """
    transcript2phylosofs = {}
    for _, gene in data.groupby('Gene stable ID'):
        for transcript_id, transcript in gene.groupby('Transcript stable ID'):
            phylosofs = []
            for subexon in transcript.sort_values(
                    'Subexon rank in transcript')['HomologousExons']:
                for ortholog in subexon.split('-'):
                    phylosofs.append(exon2char[ortholog])
            transcript2phylosofs[transcript_id] = ''.join(phylosofs)
    return transcript2phylosofs


def _get_terminal_names(input_tree):
    """Return the set of terminal names in the Newick tree at input_tree."""
    tree = Phylo.read(input_tree, "newick")
    return {clade.name for clade in tree.get_terminals()}


def _get_protein2gene(exontable_file):
    """Return a pandas datafarme to map protein (translation) id to gene id."""
    data = pd.read_csv(exontable_file, sep='\t')
    data = data.loc[:, ['Gene stable ID', 'Protein stable ID']].dropna(
    ).drop_duplicates().set_index('Protein stable ID')
    data['Gene stable ID'] = data['Gene stable ID'].astype('category')
    return data


def _get_terminals_to_delete(tree, used_proteins):
    """Return the list of terminal nodes to delete."""
    return [
        clade.name for clade in tree.get_terminals()
        if clade.name not in used_proteins and clade.name is not None
    ]


def prune_tree(input_tree, output_tree, exontable_file, used_genes):
    """
    Delete unused proteins from the Ensembl gene tree.

    The function reads the Newick tree from input_tree and stores the pruned
    tree in output_tree and returning the output_tree path. If an
    error occurs, the function does not save the tree and raises an Exception.

    In the tree, each protein is identified by their ensembl translation id.
    For each gene, Ensembl uses the longest translated transcript sequence.

    The exontable_file downloaded from Ensembl with transcript_query is needed
    for getting the mapping between translation and gene ids of the used_genes.
    """
    tree = Phylo.read(input_tree, "newick")
    data = _get_protein2gene(exontable_file)

    # Delete unused proteins
    used_proteins = set(data.index[data['Gene stable ID'].isin(used_genes)])
    for clade in _get_terminals_to_delete(tree, used_proteins):
        try:
            tree.prune(clade)
        except ValueError as err:
            raise Exception(
                'Error with protein {} while prunning the tree {} : {}'.format(
                    clade, tree, err))

    # Use gene id instead of protein/peptide/translation id in the tree
    for clade in tree.get_terminals():
        clade.name = data.loc[clade.name]

    Phylo.write(tree, output_tree, 'newick')

    return output_tree
