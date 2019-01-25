import os
import pytest
from exonhomology import transcript_info
from exonhomology import subexons


@pytest.fixture(scope='module')
def clustered_trx_data(request):
    filename = request.module.__file__

    def _get_clustered_trx_data(folder):
        test_dir = os.path.dirname(filename)
        datadir = os.path.join(test_dir, 'data')
        folder_path = os.path.join(datadir, folder)
        trx_data = transcript_info.read_transcript_info(
            os.path.join(folder_path, 'TSL', folder + '_TSL.csv'),
            os.path.join(folder_path, 'TablesExons',
                         folder + '_exonstable.tsv'),
            os.path.join(folder_path, 'Sequences', folder + '.fasta'),
            remove_na=False)
        return transcript_info.exon_clustering(trx_data)

    return _get_clustered_trx_data


def test_subexon_table(clustered_trx_data):
    data = clustered_trx_data('GRIN1_ENSG00000169258')
    subexon_table = subexons.create_subexon_table(data)

    assert 'Subexon ID' in subexon_table.columns
    assert 'Subexon ID cluster' in subexon_table.columns

    assert 'ENSMUST00000099506/ENSMUST00000135343' in subexon_table[
        'Transcript stable ID cluster'].values  # .values, otherwise in
    # checks whether the value is in the Series' index.

    assert all(subexon_table['Transcript fraction'] == 1.0)


def test_subexon_clusters(clustered_trx_data):
    data = clustered_trx_data('MAPK8_ENSG00000107643')
    subexon_table = subexons.create_subexon_table(data)
    assert len(
        subexon_table.loc[subexon_table['Subexon ID cluster'] ==
                          'ENSMMUE00000040028_SE_0/ENSMMUE00000387789_SE_0',
                          'Subexon ID'].unique()) == 1

    # 422 rows with before merging non-redundant contiguous subexons
    assert subexon_table.shape[0] == 397

    # QVQQ
    assert len(subexon_table[subexon_table['Subexon ID'] ==
                             'ENSMUSE00000689835_SE_2']) == 0
    assert len(subexon_table[subexon_table['Subexon ID'] ==
                             'ENSMUSE00000689835_SE_1']) == 0
    qvqq = subexon_table[subexon_table['Subexon ID'] ==
                         'ENSMUSE00000689835_SE_2_1']
    assert len(qvqq) == 2
    assert sorted(qvqq['Transcript stable ID']) == [
        'ENSMUST00000111943', 'ENSMUST00000111945'
    ]
    for index in [0, 1]:
        assert qvqq['Strand'].iloc[index] == -1
        assert str(qvqq['Exon protein sequence'].iloc[index]) == 'QVQQ*'
        assert str(qvqq['Subexon sequence'].iloc[index]) == 'CACAGGTGCAGCAATGA'
        assert qvqq['Exon stable ID cluster'].iloc[
            index] == 'ENSMUSE00000689835/ENSMUSE00000689841'

    # VINGSQHPSSSSSVNDVSSMSTDPTLASDTDSSLEASAGPLGCCR
    exon_row = subexon_table[subexon_table['Exon stable ID'] ==
                             'ENSMMUE00000040064'].iloc[0]
    assert exon_row['Subexon ID'] == 'ENSMMUE00000040064_SE_15_16'
    assert exon_row['Subexon ID cluster'] == 'ENSMMUE00000040064_SE_15_16'
    assert str(exon_row['Exon protein sequence']
               ) == 'VINGSQHPSSSSSVNDVSSMSTDPTLASDTDSSLEASAGPLGCCR*'
    assert str(exon_row['Subexon sequence']) == (
        'CAGTGATCAATGGCTCTCAGCATCCATCGTCATCGTCGTCTGTCAATGATGTGTCTTCAA'
        'TGTCAACAGATCCGACTTTGGCCTCGGATACAGACAGCAGTCTAGAAGCATCAGCTGGGC'
        'CTCTGGGCTGCTGTAGATGA')
    assert exon_row['Strand'] == 1
    assert exon_row['Subexon genomic coding start'] == 42539779
    assert exon_row['Subexon genomic coding end'] == 42539918
