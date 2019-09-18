import os
import pytest
from thoraxe import transcript_info
from thoraxe import subexons


@pytest.fixture(scope='module')
def subexon_data(request):
    filename = request.module.__file__

    def _get_subexon_data(folder):
        test_dir = os.path.dirname(filename)
        datadir = os.path.join(test_dir, 'data')
        folder_path = os.path.join(datadir, folder, 'Ensembl')
        trx_data = transcript_info.read_transcript_info(
            os.path.join(folder_path, 'tsl.csv'),
            os.path.join(folder_path, 'exonstable.tsv'),
            os.path.join(folder_path, 'sequences.fasta'),
            remove_na=False)
        clustered = transcript_info.exon_clustering(trx_data)
        return subexons.create_subexon_table(clustered)

    return _get_subexon_data


def test_subexon_connectivity(subexon_data):
    data = subexon_data('GPRIN1')
    connected_subexons = subexons.alignment.subexon_connectivity(data)
    # Gene: ENSMODG00000005102
    # Transcript: ENSMODT00000006412.2
    assert connected_subexons == {('ENSMODE00000054313_SE_1',
                                   'ENSMODE00000054318_SE_0')}


def test_subexon_table(subexon_data):
    data = subexon_data('GPRIN1')
    connected_subexons = subexons.alignment.subexon_connectivity(data)
    subgroup, subexon_matrix = subexons.alignment.create_subexon_matrix(data)
    chimerics = subexons.alignment.create_chimeric_sequences(
        subgroup, subexon_matrix, connected_subexons)

    assert chimerics is not None  # TO DO
