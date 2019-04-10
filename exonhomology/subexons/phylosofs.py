"""
phylosofs: Functions to generate PhyloSofS input.
"""

import string

CHARS = [
    char for char in string.printable
    if char not in {' ', '\t', '\n', '\r', '\x0b', '\x0c', '\\'}
]

MAX_EXONS = len(CHARS)


def get_exon2char(exons):
    """
    Return a dictionary from orthologous exon group to a single character.

    >>> result = exon2char(['1_0-1_1', '2_0'])
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
