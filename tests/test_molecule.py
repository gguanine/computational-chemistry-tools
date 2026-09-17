import numpy as np
import pytest

from cc_tools.molecule import Molecule


def test_create_molecule():
    mol = Molecule(
        symbols=["O", "H", "H"],
        coordinates=[
            [0.0, 0.0, 0.0],
            [0.0, 0.7, 0.5],
            [0.0, -0.7, 0.5],
        ],
        charge=0,
        multiplicity=1,
        name="water",
    )

    assert mol.symbols == ["O", "H", "H"]
    assert mol.coordinates.shape == (3, 3)
    assert mol.charge == 0
    assert mol.multiplicity == 1
    assert mol.name == "water"


def test_coordinates_are_converted_to_numpy_array():
    mol = Molecule(
        symbols=["H"],
        coordinates=[[0, 0, 0]],
    )

    assert isinstance(mol.coordinates, np.ndarray)
    assert mol.coordinates.dtype == float


def test_invalid_coordinates_raise_error():
    with pytest.raises(ValueError):
        Molecule(
            symbols=["O", "H"],
            coordinates=[[0.0, 0.0, 0.0]],
        )