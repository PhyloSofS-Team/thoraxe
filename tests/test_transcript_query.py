import glob
import os
import pytest
import shutil
import subprocess


@pytest.fixture(scope='module')
def example_folder():
    folders = glob.glob('MAPK8*')
    if len(folders) == 1:
        folder = folders[0]
        yield folder
        shutil.rmtree(folder)


def test_download():
    assert subprocess.call(["transcript_query", "MAPK8"]) == 0


def test_first_level(example_folder):
    content = glob.glob(os.path.join(example_folder, '*'))
    assert len(content) == 4
    assert any(elem.endswith('.tree.nh') for elem in content)
    assert any(elem.endswith('TSL') for elem in content)
    assert any(elem.endswith('Sequences') for elem in content)
    assert any(elem.endswith('TablesExons') for elem in content)


def test_second_level(example_folder):
    content = glob.glob(os.path.join(example_folder, '*', '*'))
    assert len(content) == 3
    assert any(elem.endswith('_TSL.csv') for elem in content)
    assert any(elem.endswith('.fasta') for elem in content)
    assert any(elem.endswith('_exonstable.tsv') for elem in content)
