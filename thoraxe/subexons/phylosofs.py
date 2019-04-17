"""
phylosofs: Functions to generate PhyloSofS input.
"""

import logging
import string

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


def prune_tree(input_path, output_path, used_proteins):
    """
    Delete unused proteins from the Ensembl gene tree.

    The function reads the Newick tree from input_path and stores the pruned
    tree in output_path if any error occurs, returning output_path. If an
    error occurs, the function does not save the tree and returns None.

    In the tree, each protein is identified by their ensembl translation id.
    For each gene, Ensembl uses the longest translated transcript sequence.
    """
    tree = Phylo.read(input_path, "newick")

    for clade in tree.get_terminals():
        if clade.name not in used_proteins:
            try:
                tree.prune(clade.name)
            except ValueError as err:
                logging.warning(
                    'ValueError with protein %s while prunning the tree: %s',
                    clade.name, err)
                return None

    Phylo.write(tree, output_path, 'newick')

    return output_path


def _get_terminal_names(input_path):
    """Return the set of terminal names in the Newick tree at input_path."""
    tree = Phylo.read(input_path, "newick")
    return {clade.name for clade in tree.get_terminals()}


# noqa TO DO : Look Gene stable ID to selected Protein stable ID before clean up data
# noqa data.loc[:,['Gene stable ID', 'Protein stable ID']
# noqa              ].drop_duplicates().groupby('Gene stable ID').apply(
# noqa      lambda x: sum(i in terminal_names for i in x['Protein stable ID']))
