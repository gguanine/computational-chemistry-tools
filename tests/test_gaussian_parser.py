from pathlib import Path
import pytest
import numpy as np

from cc_tools.gaussian import *
from cc_tools.molecule import Molecule


DATA_DIR = Path(__file__).parent / "data" / "gaussian"

def make_molecule() -> Molecule:
    return Molecule(
        symbols=["H", "H"],
        coordinates=np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.74],
        ]),
        charge=0,
        multiplicity=1,
        name="h2",
    )

def test_parse_radical_normal_gaussian_log():
    path = DATA_DIR / "radical_opt_freq_normal.log"

    mol = read_output(path).molecule

    assert mol is not None
    assert len(mol.symbols) > 0
    assert mol.coordinates.shape == (len(mol.symbols), 3)

    assert mol.charge == 0
    assert mol.multiplicity == 2
    
def test_validate_opt_freq_minimum():
    molecule = make_molecule()

    result = CalculationResult(
        molecule=molecule,
        success=True,
        frequencies=np.array([
            100.0,
            200.0,
            300.0,
        ]),
    )

    validate_opt_freq(
        result,
        require_minimum=True,
    )


def test_validate_opt_freq_rejects_imaginary():
    molecule = make_molecule()

    result = CalculationResult(
        molecule=molecule,
        success=True,
        frequencies=np.array([
            -123.4,
            100.0,
            200.0,
        ]),
    )

    with pytest.raises(
        RuntimeError,
        match="Imaginary frequencies",
    ):
        validate_opt_freq(
            result,
            require_minimum=True,
        )


def test_prepare_single_point():
    molecule = make_molecule()

    result = CalculationResult(
        molecule=molecule,
        success=True,
        frequencies=np.array([100.0]),
    )

    gjf = prepare_single_point(
        result,
        route="#p B3LYP/6-31G*",
        require_minimum=True,
        nprocs=8,
        memory="8GB",
    )

    assert gjf.molecule is molecule
    assert gjf.route == "#p B3LYP/6-31G*"
    assert gjf.nprocs == 8
    assert gjf.memory == "8GB"

def test_prepare_single_points(
    tmp_path,
    monkeypatch,
):
    input_dir = tmp_path / "opt"
    output_dir = tmp_path / "energy"

    input_dir.mkdir()

    (input_dir / "a.log").write_text("")
    (input_dir / "b.log").write_text("")

    molecule = Molecule(
        symbols=["H"],
        coordinates=[[0.0, 0.0, 0.0]],
    )

    result = CalculationResult(
        molecule=molecule,
        success=True,
        frequencies=np.array([100.0]),
    )

    monkeypatch.setattr(
        "cc_tools.gaussian.read_output",
        lambda path: result,
    )

    paths = prepare_single_points(
        input_dir,
        output_dir,
        route="#p HF/3-21G",
    )

    assert paths == [
        output_dir / "a.gjf",
        output_dir / "b.gjf",
    ]

    assert (output_dir / "a.gjf").exists()
    assert (output_dir / "b.gjf").exists()