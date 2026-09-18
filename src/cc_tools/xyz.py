# src/cc_tools/xyz.py

from __future__ import annotations

from os import PathLike
from pathlib import Path

import numpy as np

from .molecule import Molecule


def molecule_to_xyz_block(
    molecule: Molecule,
    *,
    comment: str | None = None,
) -> str:
    """Convert a Molecule to a standard XYZ block."""
    if comment is None:
        comment = molecule.name or ""

    lines = [
        str(len(molecule.symbols)),
        comment,
    ]

    for symbol, (x, y, z) in zip(
        molecule.symbols,
        molecule.coordinates,
    ):
        lines.append(
            f"{symbol:<2} {x:16.8f} {y:16.8f} {z:16.8f}"
        )

    return "\n".join(lines) + "\n"


def molecule_from_xyz_block(
    block: str,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    name: str | None = None,
    smiles: str | None = None,
) -> Molecule:
    """Create a Molecule from a standard XYZ block."""
    lines = block.splitlines()

    if len(lines) < 2:
        raise ValueError("Invalid XYZ block: expected at least two lines")

    # First line: number of atoms
    try:
        natoms = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(
            f"Invalid XYZ atom count: {lines[0]!r}"
        ) from exc

    if natoms < 1:
        raise ValueError(
            f"Invalid XYZ atom count: {natoms}"
        )

    # Second line is the XYZ comment line.
    atom_lines = lines[2:]

    if len(atom_lines) != natoms:
        raise ValueError(
            f"XYZ block declares {natoms} atoms, "
            f"but contains {len(atom_lines)} atom lines"
        )

    symbols: list[str] = []
    coordinates = np.empty((natoms, 3), dtype=float)

    for i, line in enumerate(atom_lines):
        fields = line.split()

        if len(fields) < 4:
            raise ValueError(
                f"Invalid XYZ atom line {i + 3}: {line!r}"
            )

        symbols.append(fields[0])

        try:
            coordinates[i] = [
                float(fields[1]),
                float(fields[2]),
                float(fields[3]),
            ]
        except ValueError as exc:
            raise ValueError(
                f"Invalid coordinates on XYZ line "
                f"{i + 3}: {line!r}"
            ) from exc

    return Molecule(
        symbols=symbols,
        coordinates=coordinates,
        charge=charge,
        multiplicity=multiplicity,
        name=name,
        smiles=smiles,
    )


def write_xyz(
    molecule: Molecule,
    path: str | PathLike[str],
    *,
    comment: str | None = None,
) -> None:
    """Write a Molecule to an XYZ file."""
    block = molecule_to_xyz_block(
        molecule,
        comment=comment,
    )
    Path(path).write_text(block, encoding="utf-8")


def read_xyz(
    path: str | PathLike[str],
    *,
    charge: int = 0,
    multiplicity: int = 1,
    name: str | None = None,
    smiles: str | None = None,
) -> Molecule:
    """Read a Molecule from an XYZ file."""
    block = Path(path).read_text(encoding="utf-8")

    return molecule_from_xyz_block(
        block,
        charge=charge,
        multiplicity=multiplicity,
        name=name,
        smiles=smiles,
    )