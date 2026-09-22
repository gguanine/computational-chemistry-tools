from __future__ import annotations

import shlex
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Literal
from dataclasses import dataclass
import csv

from .gaussian import GaussianInput
from .molecule import Molecule
from .xtb import optimize
from .results import CalculationResult, read_output


def prepare_gaussian_workflow(
    molecules: Sequence[Molecule],
    output_dir: str | Path,
    *,
    route_opt_freq: str,
    route_sp: str,
    use_xtb: bool = True,
    xtb_method: str = "gfn2",
    xtb_failure: Literal["error", "keep"] = "error",
    nprocs: int = 24,
    memory: str = "48GB",
    require_minimum: bool = True,
) -> Path:
    """Prepare a Gaussian Opt/Freq -> SP workflow."""

    if not molecules:
        raise ValueError("No molecules provided")

    output_dir = Path(output_dir)
    opt_dir = output_dir / "opt"
    energy_dir = output_dir / "energy"

    opt_dir.mkdir(parents=True, exist_ok=True)
    energy_dir.mkdir(parents=True, exist_ok=True)

    stems = _get_stems(molecules)

    for molecule, stem in zip(molecules, stems):
        molecule = _prepare_geometry(
            molecule,
            use_xtb=use_xtb,
            xtb_method=xtb_method,
            xtb_failure=xtb_failure,
        )

        gaussian_input = GaussianInput(
            molecule=molecule,
            route=route_opt_freq,
            nprocs=nprocs,
            memory=memory,
            checkpoint=f"{stem}.chk",
        )

        gaussian_input.write(
            opt_dir / f"{stem}.gjf"
        )

    pbs_path = output_dir / "run.pbs"

    pbs_script = render_gaussian_pbs(
        route_sp=route_sp,
        nprocs=nprocs,
        memory=memory,
        require_minimum=require_minimum,
    )

    with pbs_path.open(
        "w",
        encoding="utf-8",
        newline="\n",
    ) as f:
        f.write(pbs_script)

    return pbs_path

def _prepare_geometry(
    molecule: Molecule,
    *,
    use_xtb: bool,
    xtb_method: str,
    xtb_failure: Literal["error", "keep"],
) -> Molecule:
    if not use_xtb:
        return molecule

    try:
        return optimize(
            molecule,
            method=xtb_method,
        )

    except (RuntimeError, FileNotFoundError) as exc:
        if xtb_failure == "error":
            raise

        if xtb_failure == "keep":
            warnings.warn(
                f"xTB optimization failed for "
                f"{molecule.name or 'unnamed molecule'}; "
                f"using original geometry instead: {exc}",
                stacklevel=2,
            )
            return molecule

        raise ValueError(
            f"Unknown xTB failure policy: {xtb_failure!r}"
        )

def _get_stems(
    molecules: Sequence[Molecule],
) -> list[str]:
    stems = [
        molecule.name or f"mol_{i:04d}"
        for i, molecule in enumerate(molecules, start=1)
    ]

    if len(stems) != len(set(stems)):
        raise ValueError(
            "Molecule names must be unique within a workflow"
        )

    return stems


def render_gaussian_pbs(
    *,
    route_sp: str,
    nprocs: int = 24,
    memory: str = "48GB",
    require_minimum: bool = True,
    project: str = "gaussian",
    queue: str = "parallel",
    walltime: str = "720:00:00",
    conda_env: str = "chemistry",
) -> str:
    """Render a PBS script for Gaussian Opt/Freq -> SP calculations."""

    if nprocs < 1:
        raise ValueError("nprocs must be >= 1")

    if not route_sp.strip():
        raise ValueError("route_sp must not be empty")

    minimum_arg = (
        " --require-minimum"
        if require_minimum
        else ""
    )

    return f"""#!/bin/bash
#PBS -P {project}
#PBS -m abe
#PBS -q {queue}
#PBS -l select=1:ncpus={nprocs}:mem={memory}
#PBS -l walltime={walltime}

set -e

nprocs={nprocs}
mem={shlex.quote(memory)}

dir_opt="opt"
dir_energy="energy"

route_sp={shlex.quote(route_sp)}


# Initialize HPC environment
cd "$PBS_O_WORKDIR"

source /etc/profile.d/rec_modules.sh
module load miniconda
source ~/.bashrc
conda activate {shlex.quote(conda_env)}


run_gaussian_in_dir() {{
    local dir="$1"

    cd "$PBS_O_WORKDIR/$dir" || {{
        echo "Directory $dir not found."
        return 1
    }}

    shopt -s nullglob
    local gjf_files=(*.gjf)

    if [ ${{#gjf_files[@]}} -eq 0 ]; then
        echo "No .gjf files found in $dir"
        return 0
    fi

    for gjf_file in "${{gjf_files[@]}}"; do
        local base_name="${{gjf_file%.gjf}}"

        echo "Running Gaussian for $gjf_file"

        g16 -p="$nprocs" -m="$mem" \\
            < "$gjf_file" \\
            > "$base_name.log"

        if [ -f "$base_name.chk" ]; then
            formchk -3 \\
                "$base_name.chk" \\
                "$base_name.fchk"
        fi
    done

    cd "$PBS_O_WORKDIR"
}}


echo "=== Running Opt/Freq calculations ==="
run_gaussian_in_dir "$dir_opt"


echo "=== Preparing single-point inputs ==="
cc-tools gaussian prepare-sp \\
    "$dir_opt" \\
    "$dir_energy" \\
    --route "$route_sp" \\
    --nprocs "$nprocs" \\
    --memory "$mem" \\
    --checkpoint{minimum_arg}


echo "=== Running single-point calculations ==="
run_gaussian_in_dir "$dir_energy"


echo "=== Workflow completed ==="
"""


@dataclass
class WorkflowResult:
    """Combined result from Opt/Freq and single-point calculations.

    The Opt/Freq calculation provides thermochemical corrections.
    The single-point calculation provides the final electronic energy.

    All energies are in Hartree.
    """

    opt_freq: CalculationResult
    single_point: CalculationResult

    def __post_init__(self) -> None:
        if not self.opt_freq.success:
            raise ValueError(
                "Opt/Freq calculation was not successful"
            )

        if not self.single_point.success:
            raise ValueError(
                "Single-point calculation was not successful"
            )

        if self.opt_freq.electronic_energy is None:
            raise ValueError(
                "Opt/Freq result contains no electronic energy"
            )

        if self.single_point.electronic_energy is None:
            raise ValueError(
                "Single-point result contains no electronic energy"
            )

    @property
    def molecule(self) -> Molecule:
        """Return the optimized molecule."""
        return self.opt_freq.molecule

    @property
    def electronic_energy(self) -> float:
        """Electronic energy from the high-level single-point calculation."""
        energy = self.single_point.electronic_energy

        assert energy is not None
        return energy

    @property
    def zero_point_correction(self) -> float | None:
        """Zero-point correction from the Opt/Freq calculation."""
        return self.opt_freq.zero_point_correction

    @property
    def enthalpy_correction(self) -> float | None:
        """Thermal enthalpy correction from the Opt/Freq calculation."""
        return self.opt_freq.enthalpy_correction

    @property
    def free_energy_correction(self) -> float | None:
        """Thermal Gibbs free-energy correction from the Opt/Freq calculation."""
        return self.opt_freq.free_energy_correction

    @property
    def zero_point_energy(self) -> float | None:
        """Composite zero-point corrected energy."""
        correction = self.zero_point_correction

        if correction is None:
            return None

        return self.electronic_energy + correction

    @property
    def enthalpy(self) -> float | None:
        """Composite enthalpy."""
        correction = self.enthalpy_correction

        if correction is None:
            return None

        return self.electronic_energy + correction

    @property
    def free_energy(self) -> float | None:
        """Composite Gibbs free energy."""
        correction = self.free_energy_correction

        if correction is None:
            return None

        return self.electronic_energy + correction

    @property
    def frequencies(self):
        """Vibrational frequencies from the Opt/Freq calculation."""
        return self.opt_freq.frequencies

    @property
    def is_minimum(self) -> bool | None:
        """Whether the optimized structure is a local minimum."""
        return self.opt_freq.is_minimum

    @property
    def temperature(self) -> float | None:
        return self.opt_freq.temperature

    @property
    def pressure(self) -> float | None:
        return self.opt_freq.pressure


def collect_workflow_results(
    opt_freq_dir: str | Path,
    single_point_dir: str | Path,
    *,
    pattern: str = "*.log",
) -> list[WorkflowResult]:
    """Collect paired Opt/Freq and single-point results.

    Files are paired by identical filenames.
    """
    opt_freq_dir = Path(opt_freq_dir)
    single_point_dir = Path(single_point_dir)

    results: list[WorkflowResult] = []

    for opt_freq_path in sorted(opt_freq_dir.glob(pattern)):
        single_point_path = (
            single_point_dir / opt_freq_path.name
        )

        if not single_point_path.exists():
            raise FileNotFoundError(
                "Matching single-point output not found: "
                f"{single_point_path}"
            )

        opt_freq = read_output(opt_freq_path)
        single_point = read_output(single_point_path)

        results.append(
            WorkflowResult(
                opt_freq=opt_freq,
                single_point=single_point,
            )
        )

    return results


def workflow_result_to_dict(
    result: WorkflowResult,
) -> dict[str, object]:
    """Convert a WorkflowResult to a flat dictionary."""
    molecule = result.molecule

    imaginary_frequencies = (
        result.opt_freq.imaginary_frequencies
    )

    if imaginary_frequencies is None:
        n_imaginary = None
    else:
        n_imaginary = len(imaginary_frequencies)

    return {
        "name": molecule.name,
        "smiles": molecule.smiles,
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,

        # Calculation status
        "opt_freq_success": result.opt_freq.success,
        "optimization_converged": (
            result.opt_freq.optimization_converged
        ),
        "single_point_success": result.single_point.success,
        "is_minimum": result.is_minimum,
        "n_imaginary": n_imaginary,

        # Raw electronic energies
        "opt_freq_electronic_energy_hartree": (
            result.opt_freq.electronic_energy
        ),
        "single_point_electronic_energy_hartree": (
            result.single_point.electronic_energy
        ),

        # Thermochemical corrections
        "zero_point_correction_hartree": (
            result.zero_point_correction
        ),
        "enthalpy_correction_hartree": (
            result.enthalpy_correction
        ),
        "free_energy_correction_hartree": (
            result.free_energy_correction
        ),

        # Composite results
        "zero_point_energy_hartree": (
            result.zero_point_energy
        ),
        "enthalpy_hartree": result.enthalpy,
        "free_energy_hartree": result.free_energy,

        # Thermochemistry conditions
        "temperature_kelvin": result.temperature,
        "pressure_atm": result.pressure,
    }


def write_workflow_csv(
    results: Sequence[WorkflowResult],
    path: str | Path,
) -> None:
    """Write workflow results to a CSV file."""
    path = Path(path)

    rows = [
        workflow_result_to_dict(result)
        for result in results
    ]

    if not rows:
        raise ValueError(
            "No workflow results to write"
        )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0]),
        )

        writer.writeheader()
        writer.writerows(rows)