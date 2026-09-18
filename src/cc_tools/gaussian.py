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

@dataclass
class CalculationResult:
    molecule: Molecule
    success: bool
    frequencies: np.ndarray | None = None

def read_log(path: str | Path) -> CalculationResult:
    """Read a Gaussian log file into a CalculationResult."""
    path = Path(path)

    data = cclib.io.ccopen(str(path)).parse()

    symbols = [
        periodictable.elements[number].symbol
        for number in data.atomnos
    ]

    coordinates = np.asarray(
        data.atomcoords[-1],
        dtype=float,
    )

    molecule = Molecule(
        symbols=symbols,
        coordinates=coordinates,
        charge=int(data.charge),
        multiplicity=int(data.mult),
        name=path.stem,
    )

    frequencies = getattr(data, "vibfreqs", None)

    if frequencies is not None:
        frequencies = np.asarray(
            frequencies,
            dtype=float,
        )

    return CalculationResult(
        molecule=molecule,
        success=bool(
            data.metadata.get("success", False)
        ),
        frequencies=frequencies,
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

def validate_opt_freq(
    result: CalculationResult,
    *,
    require_minimum: bool = False,
) -> None:
    """Validate a Gaussian Opt/Freq calculation result."""
    if not result.success:
        raise RuntimeError(
            "Gaussian calculation did not terminate normally"
        )

    if result.frequencies is None:
        raise RuntimeError(
            "Gaussian calculation contains no frequencies"
        )

    if (
        require_minimum
        and np.any(result.frequencies < 0)
    ):
        imaginary = result.frequencies[
            result.frequencies < 0
        ]

        raise RuntimeError(
            "Imaginary frequencies found: "
            + ", ".join(f"{freq:.2f}" for freq in imaginary)
        )
    

def prepare_single_point(
    result: CalculationResult,
    *,
    route: str,
    require_minimum: bool = False,
    nprocs: int | None = None,
    memory: str | None = None,
    checkpoint: str | None = None,
    additional_input: str = "",
) -> GaussianInput:
    """Prepare a Gaussian single-point input from an Opt/Freq result."""
    validate_opt_freq(
        result,
        require_minimum=require_minimum,
    )

    return GaussianInput(
        molecule=result.molecule,
        route=route,
        nprocs=nprocs,
        memory=memory,
        checkpoint=checkpoint,
        additional_input=additional_input,
    )

def prepare_single_points(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    route: str,
    require_minimum: bool = False,
    nprocs: int | None = None,
    memory: str | None = None,
    checkpoint: bool = False,
    additional_input: str = "",
) -> list[Path]:
    """Prepare Gaussian single-point inputs from all log files in a directory."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    if not input_dir.is_dir():
        raise NotADirectoryError(
            f"Input directory does not exist: {input_dir}"
        )

    log_paths = sorted(input_dir.glob("*.log"))

    if not log_paths:
        raise FileNotFoundError(
            f"No .log files found in: {input_dir}"
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_paths: list[Path] = []

    for log_path in log_paths:
        result = read_log(log_path)

        stem = log_path.stem
        output_path = output_dir / f"{stem}.gjf"

        checkpoint_name = (
            f"{stem}.chk"
            if checkpoint
            else None
        )

        gaussian_input = prepare_single_point(
            result,
            route=route,
            require_minimum=require_minimum,
            nprocs=nprocs,
            memory=memory,
            checkpoint=checkpoint_name,
            additional_input=additional_input,
        )

        gaussian_input.write(output_path)

        output_paths.append(output_path)

    return output_paths