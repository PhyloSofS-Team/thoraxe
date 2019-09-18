import os
import shutil
import sys
import platform
import pytest
import thoraxe


def _is_windows():
    """Return True in Windows."""
    return platform.system() == 'Windows'


@pytest.fixture(scope='module')
def set_out_dir():
    out_dir = 'tmp'
    yield out_dir
    shutil.rmtree(out_dir)


@pytest.mark.skipif(_is_windows(), reason="Skipping this test on Windows.")
def test_thoraxe(monkeypatch, request, set_out_dir):
    filename = request.module.__file__
    in_dir = os.path.join(os.path.dirname(filename), 'data', 'MAPK8')
    monkeypatch.setattr(sys, 'argv', [
        'thoraxe', '-i', in_dir, '-o', 'tmp', '-y', '--plot_chimerics', '-a',
        'clustalo'
    ])
    assert not os.path.isdir('tmp')
    assert thoraxe.main() is None
    assert os.path.isdir('tmp')
