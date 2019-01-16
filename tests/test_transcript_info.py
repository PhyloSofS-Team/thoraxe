import os
import pytest
from exonhomology import transcript_info

@pytest.fixture(scope='module')
def datadir(request):
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    return os.path.join(test_dir, "data")


def test_read_transcript_info(datadir):
    mapk8 = os.path.join(datadir, "MAPK8_ENSG00000107643")
    trx_data = transcript_info.read_transcript_info(
        os.path.join(mapk8, "TSL", "MAPK8_ENSG00000107643_TSL.csv"),
        os.path.join(mapk8, "TablesExons", "MAPK8_ENSG00000107643_exonstable.tsv"),
        os.path.join(mapk8, "Sequences", "MAPK8_ENSG00000107643.fasta"),
        remove_na=False
        )
    # TO DO
