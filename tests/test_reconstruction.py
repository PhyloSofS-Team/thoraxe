"""
Regression tests for the protein sequences of transcripts, sub-exons and paths.

Where a test gene is used, the reference is each transcript's own CDS: its exon
sequences cut at its own genomic coding coordinates and translated,
independently of add_protein_seq. The test genes hold human transcripts only
(Ensembl 115) and were chosen because the earlier code gave them a wrong
protein.
"""
import numpy as np
import pandas as pd

from thoraxe import transcript_info


def test_identical_sequence_representative():
    """Identical proteins keep the Ensembl transcript, lower TSL, then ID."""
    data = pd.DataFrame({
        'GeneID': ['G'] * 5,
        'TranscriptID': ['CHS.1.1', 'T3', 'T2', 'T1', 'T4'],
        'Source': ['user', 'ensembl', 'havana', 'ensembl', 'ensembl'],
        'TSL': [1.0, 2.0, 1.0, np.nan, 5.0],
        'ExonProteinSequence': ['MA', 'MA', 'MA', 'MB', 'MB']
    })
    transcript_info.delete_identical_sequences(data)
    assert sorted(data['TranscriptID']) == ['T2', 'T4']
