"""
Regression tests for the protein sequences of transcripts, sub-exons and paths.

Where a test gene is used, the reference is each transcript's own CDS: its exon
sequences cut at its own genomic coding coordinates and translated,
independently of add_protein_seq. The test genes hold human transcripts only
(Ensembl 115) and were chosen because the earlier code gave them a wrong
protein.
"""
import os

import numpy as np
import pandas as pd
from Bio.Seq import Seq

from thoraxe import subexons
from thoraxe import transcript_info

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


def _read(gene):
    """Return read_transcript_info for a test gene folder."""
    folder = os.path.join(DATA, gene, 'Ensembl')
    return transcript_info.read_transcript_info(
        os.path.join(folder, 'tsl.csv'),
        os.path.join(folder, 'exonstable.tsv'),
        os.path.join(folder, 'sequences.fasta'),
        remove_na=False)


def _own_proteins(trx_data):
    """Return {TranscriptID: protein} translated from each transcript's CDS."""
    proteins = {}
    for transcript, rows in trx_data.groupby('TranscriptID'):
        cds = ''
        for row in rows.sort_values('ExonRank').itertuples():
            if row.Strand == -1:
                start = row.ExonRegionEnd - row.GenomicCodingEnd
                end = row.ExonRegionEnd - row.GenomicCodingStart + 1
            else:
                start = row.GenomicCodingStart - row.ExonRegionStart
                end = row.GenomicCodingEnd - row.ExonRegionStart + 1
            cds += str(row.ExonSequence)[int(start):int(end)]
        proteins[transcript] = str(Seq(cds).translate()).rstrip('*')
    return proteins


def _subexon_proteins(trx_data):
    """Return {TranscriptID: protein} joined from the sub-exon table."""
    table = subexons.create_subexon_table(
        transcript_info.exon_clustering(trx_data))
    table = table.sort_values(['TranscriptID', 'SubexonRank'])
    return {
        transcript: ''.join(
            str(seq) for seq in rows['SubexonProteinSequence']).replace(
                '*', '')
        for transcript, rows in table.groupby('TranscriptID')
    }


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


def test_shared_exon_keeps_own_cds():
    """Each transcript keeps its own coding bounds in a shared exon.

    OR51F2: two transcripts share their only exon with different CDS starts.
    """
    trx_data = _read('OR51F2')
    own = _own_proteins(trx_data)
    assert len(own) == 2 and len(set(own.values())) == 2
    assert _subexon_proteins(trx_data) == own
