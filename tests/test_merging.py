import os
import pytest
import numpy as np
import thoraxe as tx


@pytest.fixture(scope='module')
def kif1b(request):
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    data_dir = os.path.join(test_dir, 'data')
    kif1b_dir = os.path.join(data_dir, 'KIF1B', 'Ensembl')
    return {
        'root': os.path.join(data_dir, 'KIF1B'),
        'tsl': os.path.join(kif1b_dir, 'tsl.csv'),
        'exontable': os.path.join(kif1b_dir, 'exonstable.tsv'),
        'seqs': os.path.join(kif1b_dir, 'sequences.fasta')
    }


def test_read_transcript_info(kif1b):
    transcript_table = tx.get_transcripts(kif1b['root'])

    # there is only data from human and mouse
    assert all(
        transcript_table.Species.unique() == ['homo_sapiens', 'mus_musculus'])

    exon_a = 'ENSE00000955614'
    exon_b = 'ENSE00000955628'

    rows_a = transcript_table.ExonIDCluster == exon_a
    rows_b = transcript_table.ExonIDCluster == exon_b

    assert sorted(transcript_table.loc[rows_a,
                                       'TranscriptIDCluster'].unique()) == [
                                           'ENST00000263934',
                                           'ENST00000377083/ENST00000377093',
                                           'ENST00000622724'
                                       ]
    assert sorted(transcript_table.loc[rows_b,
                                       'TranscriptIDCluster'].unique()) == [
                                           'ENST00000377081',
                                           'ENST00000377086', 'ENST00000620295'
                                       ]

    # the two exons share the 3' region:

    assert all(
        transcript_table.loc[rows_a,
                             'GenomicCodingStart'].unique() == [10282400])
    assert all(transcript_table.loc[rows_a,
                                    'GenomicCodingEnd'].unique() == [10282533])

    assert all(
        transcript_table.loc[rows_b,
                             'GenomicCodingStart'].unique() == [10282322])
    assert all(transcript_table.loc[rows_b,
                                    'GenomicCodingEnd'].unique() == [10282533])

    # there is not differences at the nucleotide level:

    assert len(transcript_table.loc[rows_a,
                                    'ExonSequence'].map(str).unique()) == 1
    assert len(transcript_table.loc[rows_b,
                                    'ExonSequence'].map(str).unique()) == 1

    # but the residue at the beginning changes depending on the transcript
    # i.e. there is a shared residue with another exon:

    assert len(
        transcript_table.loc[rows_a,
                             'ExonProteinSequence'].map(str).unique()) == 2
    assert len(
        transcript_table.loc[rows_b,
                             'ExonProteinSequence'].map(str).unique()) == 2

    clustered = tx.transcript_info.exon_clustering(
        transcript_table,
        minimum_len=4,
        coverage_cutoff=80.0,
        percent_identity_cutoff=30.0,
        gap_open_penalty=-10,
        gap_extend_penalty=-1)

    # all belong to the same cluster:

    cluster = set(clustered.loc[clustered.ExonIDCluster == exon_a, 'Cluster'])
    cluster.update(clustered.loc[clustered.ExonIDCluster == exon_b, 'Cluster'])
    assert len(cluster) == 1

    subexon_table = tx.subexons.create_subexon_table(clustered)

    # subexons 13 and 14 shouldn't be merged, because the two 14 should be merged first

    assert sorted(
        subexon_table.loc[subexon_table.ExonIDCluster == exon_a,
                          'SubexonIDCluster'].unique()) == [
                              'ENSE00000955614_SE_14/ENSE00000955628_SE_14'
                          ]
    assert sorted(
        subexon_table.loc[subexon_table.ExonIDCluster == exon_b,
                          'SubexonIDCluster'].unique()) == [
                              'ENSE00000955614_SE_14/ENSE00000955628_SE_14',
                              'ENSE00000955628_SE_13'
                          ]
