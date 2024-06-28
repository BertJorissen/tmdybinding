"""Tight-binding lattices for TMDs for the Pybinding-package."""
from .tmd_lattice import *
from .tmd_abstract_lattice import *
from .parameters import *


def tests():
    """Run the package tests."""
    import pytest
    import pathlib
    import os
    from pybinding.utils.misc import cd
    import pybinding as pb

    module_path = pathlib.Path(__file__).parent

    args = []
    if (module_path / 'tests').exists():
        # tests are inside installed package -> use read-only mode
        args.append('--failpath=' + os.getcwd() + '/failed')
        with cd(module_path):
            args += ['-c', str(module_path / 'tests/local.cfg'), str(module_path)]
            error_code = pytest.main(args)
    else:
        # tests are in dev environment -> use development mode
        with cd(module_path.parent):
            error_code = pytest.main(args)

    return error_code or None
