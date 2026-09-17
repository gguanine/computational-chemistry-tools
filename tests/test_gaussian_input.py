import numpy as np
from cc_tools.molecule import Molecule
from cc_tools.gaussian import GaussianInput, read_gjf


def test_gaussian_input_render():
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

    gjf = GaussianInput(
        molecule=mol,
        route="#p wb97xd/def2svp opt freq",
        nprocs=8,
        memory="16GB",
        checkpoint="water.chk",
    )

    text = gjf.render()

    assert "%nprocshared=8" in text
    assert "%mem=16GB" in text
    assert "%chk=water.chk" in text

    assert "#p wb97xd/def2svp opt freq" in text

    assert "water" in text
    assert "0 1" in text

    assert "O" in text
    assert "H" in text

def test_gaussian_input_without_link0_commands():
    mol = Molecule(
        symbols=["He"],
        coordinates=[[0.0, 0.0, 0.0]],
        name="helium",
    )

    gjf = GaussianInput(
        molecule=mol,
        route="#p hf/sto-3g",
    )

    text = gjf.render()

    assert "%nprocshared" not in text
    assert "%mem" not in text
    assert "%chk" not in text

    assert "#p hf/sto-3g" in text

def test_gaussian_input_write(tmp_path):
    mol = Molecule(
        symbols=["He"],
        coordinates=[[0.0, 0.0, 0.0]],
        name="helium",
    )

    gjf = GaussianInput(
        molecule=mol,
        route="#p hf/sto-3g",
    )

    output = tmp_path / "helium.gjf"

    gjf.write(output)

    assert output.exists()

    text = output.read_text()

    assert "#p hf/sto-3g" in text
    assert "helium" in text

def test_gaussian_input_round_trip(tmp_path):
    original = Molecule(
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

    path = tmp_path / "water.gjf"

    GaussianInput(
        molecule=original,
        route="#p hf/sto-3g",
    ).write(path)

    loaded = read_gjf(path)

    assert loaded.name == original.name
    assert loaded.symbols == original.symbols
    assert loaded.charge == original.charge
    assert loaded.multiplicity == original.multiplicity

    np.testing.assert_allclose(
        loaded.coordinates,
        original.coordinates,
    )