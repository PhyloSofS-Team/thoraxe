"""
split: Util functions to split homologous exons.
"""

import numpy as np


def split_exons(exons):
    """
    Return a list of the homologous exons in the string.

    >>> split_exons('3_0/3_1')
    ['3_0', '3_1']
    >>> split_exons('4_0')
    ['4_0']
    """
    if '/' in exons:
        return exons.split('/')

    return [exons]


def _int_length(length, length_seq):
    """
    Return length as an int.

    >>> _int_length("", 10.0)
    10
    """
    if not isinstance(length, str) and np.isnan(length):
        return length_seq

    if length == '':
        return int(length_seq)

    return int(length)


def split_lengths(lengths, seq):
    """
    Return a list of the homologous exon lengths, len(seq) if there is not data.

    >>> import numpy as np
    >>> split_lengths('4/2', 'MLGHAC')
    [4, 2]
    >>> split_lengths('6', 'MLGHAC')
    [6]
    >>> split_lengths('', 'MLGHAC')
    [6]
    >>> split_lengths(np.nan, 'MLGHAC')
    [6]
    """
    length_seq = len(seq)

    if isinstance(lengths, str) and '/' in lengths:
        return [int(length) for length in lengths.split('/')]

    return [_int_length(lengths, length_seq)]


def split_seqs(sequences, seq):
    """
    Return a list with the homologous exon sequences.

    >>> split_seqs('MLGH/AC', 'MLGHAC')
    ['MLGH', 'AC']
    >>> split_seqs('MLGHAC', 'MLGHAC')
    ['MLGHAC']
    >>> split_seqs('', 'MLGHAC')
    ['MLGHAC']
    """
    if '/' in sequences:
        return sequences.split('/')

    if not sequences:
        return [seq]

    return [sequences]
