import os
import shutil
import sys
import platform
import subprocess
import pytest

import thoraxe
from thoraxe import utils


@pytest.fixture(scope='module')
def set_out_dir():
    out_dir = 'tmp'
    yield out_dir
    os.path.isdir(out_dir) and shutil.rmtree(out_dir)


def test_thoraxe(monkeypatch, request, set_out_dir):
    aligner = 'ProGraphMSA'
    if utils.windows.is_windows():
        if shutil.which('wsl') is not None:  # ProGraphSA in the WSL
            status, _ = subprocess.getstatusoutput(
                ['wsl', 'ProGraphMSA', '-h'])
            if status == 0:
                aligner = 'wsl ProGraphMSA'
            else:
                raise Exception(
                    'You should install ProGraphMSA in the Windows '
                    'Subsystem for Linux to run this test.')
        else:  # AppVeyor
            aligner = ""

    if aligner:
        filename = request.module.__file__
        in_dir = os.path.join(os.path.dirname(filename), 'data', 'MAPK8')
        monkeypatch.setattr(sys, 'argv', [
            'thoraxe', '-i', in_dir, '-o', 'tmp', '-y', '--plot_chimerics',
            '-a', aligner
        ])
        assert not os.path.isdir('tmp')
        assert thoraxe.main() is None
        assert os.path.isdir('tmp')
