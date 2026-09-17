from dataclasses import dataclass
from pathlib import Path

import cclib
import numpy as np
import periodictable

from .molecule import Molecule


@dataclass
class GaussianInput:
    molecule: Molecule
    route: str

    nprocs: int | None = None
    memory: str | None = None
    checkpoint: str | None = None

    additional_input: str = ""

    def render(self) -> str:
        mol = self.molecule

        lines = []

        if self.nprocs is not None:
            lines.append(f"%nprocshared={self.nprocs}")

        if self.memory is not None:
            lines.append(f"%mem={self.memory}")

        if self.checkpoint is not None:
            lines.append(f"%chk={self.checkpoint}")

        lines.extend([
            self.route,
            "",
            mol.name or "Untitled",
            "",
            f"{mol.charge} {mol.multiplicity}",
        ])

        for symbol, xyz in zip(mol.symbols, mol.coordinates):
            x, y, z = xyz
            lines.append(
                f"{symbol:<3} {x:15.8f} {y:15.8f} {z:15.8f}"
            )

        lines.append("")

        if self.additional_input:
            lines.append(self.additional_input.rstrip())
            lines.append("")

        return "\n".join(lines) + "\n"

    def write(self, path: str | Path) -> None:
        Path(path).write_text(self.render())

def molecule_from_log(
    path: str | Path,
    *,
    name_suffix: str = "out",
    allow_imaginary: bool = True,
) -> Molecule:
    path = Path(path)

    data = cclib.io.ccopen(str(path)).parse()

    if not data.metadata.get("success", False):
        raise RuntimeError(
            f"Gaussian did not terminate normally: {path}"
        )

    frequencies = getattr(data, "vibfreqs", None)

    if (
        frequencies is not None
        and len(frequencies)
        and frequencies[0] < 0
        and not allow_imaginary
    ):
        raise RuntimeError(
            f"Imaginary frequency found: {path}"
        )

    symbols = [
        periodictable.elements[number].symbol
        for number in data.atomnos
    ]

    coordinates = np.asarray(
        data.atomcoords[-1],
        dtype=float,
    )

    return Molecule(
        symbols=symbols,
        coordinates=coordinates,
        charge=int(data.charge),
        multiplicity=int(data.mult),
        name=f"{path.stem}_{name_suffix}",
    )

def read_gjf(path: str | Path) -> Molecule:
    text = Path(path).read_text().strip()

    sections = text.split("\n\n")

    if len(sections) < 3:
        raise ValueError(f"Invalid Gaussian input: {path}")

    name = sections[1].strip() or None

    geometry_lines = sections[2].splitlines()

    charge, multiplicity = map(
        int,
        geometry_lines[0].split()[:2],
    )

    symbols = []
    coordinates = []

    for line in geometry_lines[1:]:
        parts = line.split()

        if len(parts) < 4:
            continue

        symbols.append(parts[0])
        coordinates.append(
            [float(x) for x in parts[1:4]]
        )

    return Molecule(
        symbols=symbols,
        coordinates=coordinates,
        charge=charge,
        multiplicity=multiplicity,
        name=name,
    )