import os
import pytest

from thoraxe import transcript_info
from thoraxe import subexons


@pytest.fixture(scope='module')
def clustered_trx_data(request):
    filename = request.module.__file__

    def _get_clustered_trx_data(folder):
        test_dir = os.path.dirname(filename)
        datadir = os.path.join(test_dir, 'data')
        folder_path = os.path.join(datadir, folder, 'Ensembl')
        trx_data = transcript_info.read_transcript_info(
            os.path.join(folder_path, 'tsl.csv'),
            os.path.join(folder_path, 'exonstable.tsv'),
            os.path.join(folder_path, 'sequences.fasta'),
            remove_na=False)
        return transcript_info.exon_clustering(trx_data)

    return _get_clustered_trx_data


def test_subexon_table(clustered_trx_data):
    data = clustered_trx_data('GPRIN1')
    subexon_table = subexons.create_subexon_table(data)

    assert 'SubexonID' in subexon_table.columns
    assert 'SubexonIDCluster' in subexon_table.columns

    assert 'ENSMUST00000099506/ENSMUST00000135343' in subexon_table[
        'TranscriptIDCluster'].values  # .values, otherwise in
    # checks whether the value is in the Series' index.

    assert all(subexon_table['TranscriptFraction'] == 1.0)


def test_subexon_clusters(clustered_trx_data):
    data = clustered_trx_data('MAPK8')
    subexon_table = subexons.create_subexon_table(data)
    assert len(
        subexon_table.loc[subexon_table['SubexonIDCluster'] ==
                          'ENSMMUE00000040028_SE_0/ENSMMUE00000387789_SE_0',
                          'SubexonID'].unique()) == 1

    not_merged = subexons.create_subexon_table(data, merge_non_redundant=False)
    assert not_merged.shape[0] > subexon_table.shape[0]

    # QVQQ
    assert len(
        not_merged[not_merged['SubexonID'] == 'ENSMUSE00000689835_SE_2']) > 0
    assert len(
        not_merged[not_merged['SubexonID'] == 'ENSMUSE00000689835_SE_1']) > 0
    assert len(subexon_table[subexon_table['SubexonID'] ==
                             'ENSMUSE00000689835_SE_2']) == 0
    assert len(subexon_table[subexon_table['SubexonID'] ==
                             'ENSMUSE00000689835_SE_1']) == 0
    qvqq = subexon_table[subexon_table['SubexonID'] ==
                         'ENSMUSE00000689835_SE_2_1']
    assert len(qvqq) == 2
    assert sorted(
        qvqq['TranscriptID']) == ['ENSMUST00000111943', 'ENSMUST00000111945']
    for index in [0, 1]:
        assert qvqq['Strand'].iloc[index] == -1
        assert str(qvqq['SubexonProteinSequence'].iloc[index]) == 'QVQQ*'
        assert str(qvqq['SubexonSequence'].iloc[index]) == 'CACAGGTGCAGCAATGA'
        assert qvqq['ExonIDCluster'].iloc[
            index] == 'ENSMUSE00000689835/ENSMUSE00000689841'

    # VINGSQHPSSSSSVNDVSSMSTDPTLASDTDSSLEASAGPLGCCR
    exon_row = subexon_table[subexon_table['ExonID'] ==
                             'ENSMMUE00000040064'].iloc[0]
    assert exon_row['SubexonID'] == 'ENSMMUE00000040064_SE_15_16'
    assert exon_row['SubexonIDCluster'] == 'ENSMMUE00000040064_SE_15_16'
    assert str(exon_row['SubexonProteinSequence']
               ) == 'AVINGSQHPSSSSSVNDVSSMSTDPTLASDTDSSLEASAGPLGCCR*'
    assert str(exon_row['SubexonSequence']) == (
        'CAGTGATCAATGGCTCTCAGCATCCATCGTCATCGTCGTCTGTCAATGATGTGTCTTCAA'
        'TGTCAACAGATCCGACTTTGGCCTCGGATACAGACAGCAGTCTAGAAGCATCAGCTGGGC'
        'CTCTGGGCTGCTGTAGATGA')
    assert exon_row['Strand'] == 1
    assert exon_row['SubexonCodingStart'] == 42539779
    assert exon_row['SubexonCodingEnd'] == 42539918
