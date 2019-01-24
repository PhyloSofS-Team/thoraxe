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
