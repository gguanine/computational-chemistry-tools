# tests/integration/test_xtb.py

import shutil

import numpy as np
import pytest

from cc_tools.molecule import Molecule
from cc_tools.xtb import optimize


pytestmark = pytest.mark.integration


@pytest.mark.skipif(
    shutil.which("xtb") is None,
    reason="xTB executable not available",
)
def test_xtb_optimize_real():
    mol = Molecule(
        symbols=["O", "H", "H"],
        coordinates=np.array([
            [0.000, 0.000, 0.000],
            [0.750, 0.000, 0.500],
            [-0.750, 0.000, 0.500],
        ]),
        charge=0,
        multiplicity=1,
        name="water",
        smiles="O",
    )

    result = optimize(mol)

    assert result.symbols == mol.symbols
    assert result.coordinates.shape == (3, 3)
    assert np.isfinite(result.coordinates).all()

    # Optimization should have changed the geometry.
    assert not np.allclose(
        result.coordinates,
        mol.coordinates,
    )