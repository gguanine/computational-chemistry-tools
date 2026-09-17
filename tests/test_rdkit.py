import numpy as np
import pytest

from cc_tools.rdkit import from_smiles


def test_from_smiles():
    mol = from_smiles(
        "O",
        name="water",
    )

    assert mol.name == "water"
    assert mol.smiles == "O"

    assert mol.symbols.count("O") == 1
    assert mol.symbols.count("H") == 2

    assert mol.charge == 0
    assert mol.multiplicity == 1

    assert mol.coordinates.shape == (3, 3)


def test_from_smiles_generates_3d_coordinates():
    mol = from_smiles("CCO")

    assert isinstance(mol.coordinates, np.ndarray)
    assert mol.coordinates.shape[1] == 3


def test_invalid_smiles():
    with pytest.raises(ValueError):
        from_smiles("this-is-not-smiles")