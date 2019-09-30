import glob
import os
import shutil
import sys
import pytest

from thoraxe import transcript_query


def _is_travis():
    """Return True in Travis CI."""
    return "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


@pytest.fixture(scope='module')
def example_folder():
    folders = glob.glob('MAPK8')
    if len(folders) == 1:
        folder = folders[0]
        yield folder
        shutil.rmtree(folder)


@pytest.mark.skipif(_is_travis(), reason="Skipping this test on Travis CI.")
def test_download(monkeypatch):
    monkeypatch.setattr(
        sys,
        'argv',
        [
            'transcript_query',
            'MAPK8',
            '-l',
            'homo_sapiens,mus_musculus,pan_troglodytes,'
            'panthera_tigris_altaica,cebus_capucinus_imitator,'
            # 'cricetulus_griseus_chok1gshd,'
            # 'colobus_angolensis_palliatus,'
            # 'cricetulus_griseus_crigri,'
            # 'heterocephalus_glaber_female,heterocephalus_glaber_male,'
            # 'canis_lupus_dingo,''
            'mustela_putorius_furo'
        ])
    assert not os.path.isdir('MAPK8')
    assert transcript_query.main() is None
    assert os.path.isdir('MAPK8')


@pytest.mark.skipif(_is_travis(), reason="Skipping this test on Travis CI.")
def test_first_level(example_folder):
    content = glob.glob(os.path.join(example_folder, '*'))
    assert len(content) == 1
    assert any(elem.endswith('Ensembl') for elem in content)


@pytest.mark.skipif(_is_travis(), reason="Skipping this test on Travis CI.")
def test_second_level(example_folder):
    content = glob.glob(os.path.join(example_folder, '*', '*'))
    assert len(content) == 4
    assert any(elem.endswith('tree.nh') for elem in content)
    assert any(elem.endswith('tsl.csv') for elem in content)
    assert any(elem.endswith('.fasta') for elem in content)
    assert any(elem.endswith('exonstable.tsv') for elem in content)

    # MAPK8/Ensembl/exonstable.tsv
    # Query ERROR: caught BioMart::Exception: non-BioMart die():
    for filename in content:
        if filename.endswith('exonstable.tsv'):
            with open(filename, 'r') as file:
                data = file.read()
                assert "ERROR" not in data
