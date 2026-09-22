from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import cclib
import numpy as np
import periodictable

from .molecule import Molecule


EV_PER_HARTREE = 27.21138505


@dataclass
class CalculationResult:
    """Parsed result of a quantum chemistry calculation.

    Energy values are stored in Hartree.
    Frequencies are stored in cm^-1.
    Temperature is stored in K.
    Pressure is stored in atm.
    """

    molecule: Molecule
    success: bool

    # Electronic energy
    electronic_energy: float | None = None

    # Geometry optimization
    optimization_converged: bool | None = None

    # Vibrational analysis
    frequencies: np.ndarray | None = None

    # Thermochemistry
    zero_point_energy: float | None = None
    enthalpy: float | None = None
    free_energy: float | None = None
    entropy: float | None = None

    temperature: float | None = None
    pressure: float | None = None

    @property
    def imaginary_frequencies(self) -> np.ndarray | None:
        """Return imaginary frequencies, if frequencies are available."""
        if self.frequencies is None:
            return None

        return self.frequencies[self.frequencies < 0]

    @property
    def is_minimum(self) -> bool | None:
        """Whether the frequency calculation corresponds to a minimum."""
        if self.frequencies is None:
            return None

        return not np.any(self.frequencies < 0)

    @property
    def zero_point_correction(self) -> float | None:
        """Return the zero-point correction relative to electronic energy."""
        if (
            self.zero_point_energy is None
            or self.electronic_energy is None
        ):
            return None

        return self.zero_point_energy - self.electronic_energy

    @property
    def enthalpy_correction(self) -> float | None:
        """Return the thermal enthalpy correction."""
        if (
            self.enthalpy is None
            or self.electronic_energy is None
        ):
            return None

        return self.enthalpy - self.electronic_energy

    @property
    def free_energy_correction(self) -> float | None:
        """Return the thermal Gibbs free-energy correction."""
        if (
            self.free_energy is None
            or self.electronic_energy is None
        ):
            return None

        return self.free_energy - self.electronic_energy


def read_output(path: str | Path) -> CalculationResult:
    """Read a quantum chemistry output file using cclib.

    The file format is detected by cclib, so this function can be used
    for supported programs such as Gaussian and ORCA.

    Electronic energies reported by cclib are converted from eV to
    Hartree. Thermochemical quantities are already reported by cclib
    in Hartree/particle.
    """
    path = Path(path)

    data = cclib.io.ccread(str(path))

    if data is None:
        raise ValueError(
            f"Could not parse quantum chemistry output: {path}"
        )

    molecule = _read_molecule(
        data,
        name=path.stem,
    )

    return CalculationResult(
        molecule=molecule,
        success=_read_success(data),
        electronic_energy=_read_electronic_energy(data),
        optimization_converged=_read_optimization_converged(data),
        frequencies=_read_array(data, "vibfreqs"),
        zero_point_energy=_read_float(data, "zpve"),
        enthalpy=_read_float(data, "enthalpy"),
        free_energy=_read_float(data, "freeenergy"),
        entropy=_read_float(data, "entropy"),
        temperature=_read_float(data, "temperature"),
        pressure=_read_float(data, "pressure"),
    )


def _read_molecule(
    data: Any,
    *,
    name: str | None = None,
) -> Molecule:
    """Create a Molecule from parsed cclib data."""
    atomnos = getattr(data, "atomnos", None)
    atomcoords = getattr(data, "atomcoords", None)

    if atomnos is None:
        raise ValueError(
            "Parsed output contains no atomic numbers"
        )

    if atomcoords is None or len(atomcoords) == 0:
        raise ValueError(
            "Parsed output contains no molecular geometry"
        )

    symbols = [
        periodictable.elements[int(number)].symbol
        for number in atomnos
    ]

    coordinates = np.asarray(
        atomcoords[-1],
        dtype=float,
    )

    charge = getattr(data, "charge", 0)
    multiplicity = getattr(data, "mult", 1)

    return Molecule(
        symbols=symbols,
        coordinates=coordinates,
        charge=int(charge),
        multiplicity=int(multiplicity),
        name=name,
    )


def _read_electronic_energy(
    data: Any,
) -> float | None:
    """Return the best available electronic energy in Hartree.

    cclib stores electronic energies in eV.

    Priority:
        CC > MP > SCF
    """
    for attribute in (
        "ccenergies",
        "mpenergies",
        "scfenergies",
    ):
        energies = getattr(data, attribute, None)

        if energies is None or len(energies) == 0:
            continue

        energy = np.asarray(energies)[-1]

        # mpenergies and ccenergies can be rank-2 arrays.
        # The last element represents the highest available
        # correction level for the final geometry.
        energy = np.asarray(energy).reshape(-1)[-1]

        return float(energy) / EV_PER_HARTREE

    return None


def _read_success(data: Any) -> bool:
    """Return whether the calculation terminated successfully."""
    metadata = getattr(data, "metadata", {})

    return bool(metadata.get("success", False))


def _read_optimization_converged(
    data: Any,
) -> bool | None:
    """Return geometry optimization convergence status."""
    optdone = getattr(data, "optdone", None)

    if optdone is None:
        return None

    return bool(optdone)


def _read_float(
    data: Any,
    attribute: str,
) -> float | None:
    """Read an optional scalar float attribute."""
    value = getattr(data, attribute, None)

    if value is None:
        return None

    return float(value)


def _read_array(
    data: Any,
    attribute: str,
) -> np.ndarray | None:
    """Read an optional NumPy array attribute."""
    value = getattr(data, attribute, None)

    if value is None:
        return None

    return np.asarray(
        value,
        dtype=float,
    )