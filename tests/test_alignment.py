import os
import pytest
from exonhomology import transcript_info
from exonhomology import subexons


@pytest.fixture(scope='module')
def subexon_data(request):
    filename = request.module.__file__

    def _get_subexon_data(folder):
        test_dir = os.path.dirname(filename)
        datadir = os.path.join(test_dir, 'data')
        folder_path = os.path.join(datadir, folder)
        trx_data = transcript_info.read_transcript_info(
            os.path.join(folder_path, 'TSL', folder + '_TSL.csv'),
            os.path.join(folder_path, 'TablesExons',
                         folder + '_exonstable.tsv'),
            os.path.join(folder_path, 'Sequences', folder + '.fasta'),
            remove_na=False)
        clustered = transcript_info.exon_clustering(trx_data)
        return subexons.create_subexon_table(clustered)

    return _get_subexon_data


def test_subexon_table(subexon_data):
    data = subexon_data('GRIN1_ENSG00000169258')
    assert True
