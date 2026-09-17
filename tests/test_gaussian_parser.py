from pathlib import Path

from cc_tools.gaussian import molecule_from_log


DATA_DIR = Path(__file__).parent / "data" / "gaussian"


def test_parse_radical_normal_gaussian_log():
    path = DATA_DIR / "radical_opt_freq_normal.log"

    mol = molecule_from_log(path)

    assert mol is not None
    assert len(mol.symbols) > 0
    assert mol.coordinates.shape == (len(mol.symbols), 3)

    assert mol.charge == 0
    assert mol.multiplicity == 2