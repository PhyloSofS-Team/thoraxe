import os
import shutil
import sys
import platform
import pytest
import subprocess
import warnings

import thoraxe


def _is_windows():
    """Return True in Windows."""
    return platform.system() == 'Windows'


@pytest.fixture(scope='module')
def set_out_dir():
    out_dir = 'tmp'
    yield out_dir
    shutil.rmtree(out_dir)


def test_thoraxe(monkeypatch, request, set_out_dir):
    aligner = 'clustalo'
    if _is_windows():
        if shutil.which('wsl') is not None:  # Clustal Omega in the WSL
            status, _ = subprocess.getstatusoutput(['wsl', 'clustalo', '-h'])
            if status == 0:
                aligner = 'wsl clustalo'
            else:
                raise Exception('You should install clustalo in the Windows '
                                'Subsystem for Linux to run this test.')
        else:  # AppVeyor
            muscle = 'C:\\projects\\thoraxe\\muscle.exe'
            warnings.warn('Assuming that the test is running in AppVeyor with '
                          'muscle installed at {}'.format(muscle))
            if os.path.exists(muscle):
                aligner = muscle + ' -in '

    filename = request.module.__file__
    in_dir = os.path.join(os.path.dirname(filename), 'data', 'MAPK8')
    monkeypatch.setattr(sys, 'argv', [
        'thoraxe', '-i', in_dir, '-o', 'tmp', '-y', '--plot_chimerics', '-a',
        aligner
    ])
    assert not os.path.isdir('tmp')
    assert thoraxe.main() is None
    assert os.path.isdir('tmp')
