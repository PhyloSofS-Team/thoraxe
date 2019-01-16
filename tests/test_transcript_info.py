import os
import pytest
from exonhomology import transcript_info

@pytest.fixture(scope='module')
def datadir(request):
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    return os.path.join(test_dir, "data")


@pytest.fixture(scope='module')
def mapk8(datadir):
    mapk8_dir os.path.join(datadir, "MAPK8_ENSG00000107643")
    return {'tsl' : os.path.join(mapk8, "TSL", "MAPK8_ENSG00000107643_TSL.csv"),},
            'exontable' : os.path.join(mapk8, "TablesExons", "MAPK8_ENSG00000107643_exonstable.tsv"),
            'seqs' : os.path.join(mapk8, "Sequences", "MAPK8_ENSG00000107643.fasta")}


def test_read_transcript_info(mapk8):
    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'],
        mapk8['exontable'],
        mapk8['seqs']
        )

    assert trx_data.loc[trx_data['Transcript stable ID'] == 'ENST00000374179',
        'Flags'].unique()[0] == '1 (assigned to previous version 7)'

    # Only h. sapiens & m. musculus have TSL information:
    assert all(trx_data.Species.unique() == ['homo_sapiens', 'mus_musculus'])


def test_remove_na(mapk8):
    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'],
        mapk8['exontable'],
        mapk8['seqs'],
        remove_na=False
        )

    # I keep other species, not only h. sapiens & m. musculus:
    assert len(trx_data.Species.unique()) > 2

    # ENSRNOT00000083933 has Xs in its sequence: ...VILGMGYKENGQXVXHVQRGLICC*
    assert 0 == sum(trx_data['Transcript stable ID'] == 'ENSRNOT00000083933')


def test_keep_badquality_sequences(mapk8):
    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'],
        mapk8['exontable'],
        mapk8['seqs'],
        remove_na=False,
        remove_badquality=False
        )

    # ENSRNOT00000083933 has Xs in its sequence: ...VILGMGYKENGQXVXHVQRGLICC*
    assert 5 == sum(trx_data['Transcript stable ID'] == 'ENSRNOT00000083933')
