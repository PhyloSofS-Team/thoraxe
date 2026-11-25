"""
Regression tests for the protein sequences of transcripts, sub-exons and paths.

Where a test gene is used, the reference is each transcript's own CDS: its exon
sequences cut at its own genomic coding coordinates and translated,
independently of add_protein_seq. The test genes hold human transcripts only
(Ensembl 115) and were chosen because the earlier code gave them a wrong
protein.
"""
import os
import shutil
import sys

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq

import thoraxe
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


def test_merged_subexons_keep_own_protein():
    """Each transcript keeps its own residues in a merged sub-exon.

    LBX2: a codon split at a merged sub-exon boundary differs between
    transcripts.
    """
    trx_data = _read('LBX2')
    assert _subexon_proteins(trx_data) == _own_proteins(trx_data)


@pytest.mark.skipif(shutil.which('ProGraphMSA') is None,
                    reason='ProGraphMSA is not installed')
def test_path_protein_is_own_protein(monkeypatch, tmp_path):
    """The s-exon sequences along each path give the transcript's protein.

    LBX2 (see test_merged_subexons_keep_own_protein).
    """
    out_dir = str(tmp_path / 'LBX2')
    monkeypatch.setattr(
        sys, 'argv',
        ['thoraxe', '-i',
         os.path.join(DATA, 'LBX2'), '-o', out_dir, '-y'])
    assert thoraxe.main() is None
    table = pd.read_csv(os.path.join(out_dir, 'thoraxe', 's_exon_table.csv'),
                        keep_default_na=False)
    table = table.sort_values(
        ['TranscriptIDCluster', 'SubexonRank', 'S_exon_Rank'], kind='stable')
    paths = {}
    for cluster, rows in table.groupby('TranscriptIDCluster'):
        protein = ''.join(rows['S_exon_Sequence']).replace('*', '')
        for transcript in cluster.split('/'):
            paths[transcript] = protein
    own = _own_proteins(_read('LBX2'))
    assert {transcript: paths[transcript] for transcript in own} == own


def test_normalize_user_phases():
    """User-added transcripts get phase -1 at their UTR boundaries."""
    data = pd.DataFrame({
        'TranscriptID': ['CHS.1.1'] * 3 + ['CHS.2.1'] + ['ENST1'] * 2,
        'Source': ['user'] * 4 + ['ensembl'] * 2,
        'Strand': [1, 1, 1, -1, 1, 1],
        'ExonRank': [1, 2, 3, 1, 1, 2],
        'ExonRegionStart': [100, 300, 500, 900, 100, 300],
        'ExonRegionEnd': [200, 400, 600, 1000, 200, 400],
        'GenomicCodingStart': [150, 300, 500, 920, 150, 300],
        'GenomicCodingEnd': [200, 400, 550, 950, 200, 350],
        'StartPhase': [0, 0, 2, 0, 0, 0],
        'EndPhase': [0, 2, 0, 0, 0, 0]
    })
    transcript_info.normalize_user_phases(data)
    assert data['StartPhase'].tolist() == [-1, 0, 2, -1, 0, 0]
    assert data['EndPhase'].tolist() == [0, 2, -1, -1, 0, 0]


def test_single_coding_exon_with_utr():
    """The only coding exon is cut at both ends of its CDS.

    OR5B21: the only coding exon has a 3'UTR (phases 0/-1).
    """
    trx_data = _read('OR5B21')
    own = _own_proteins(trx_data)
    assert list(own) == ['ENST00000360374']
    protein = ''.join(str(seq) for seq in trx_data['ExonProteinSequence'])
    assert protein.replace('*', '') == own['ENST00000360374']
