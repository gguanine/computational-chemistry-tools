from __future__ import annotations

import shlex
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

from .gaussian import GaussianInput
from .molecule import Molecule
from .xtb import optimize


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

    pbs_path.write_text(
        render_gaussian_pbs(
            route_sp=route_sp,
            nprocs=nprocs,
            memory=memory,
            require_minimum=require_minimum,
        ),
        encoding="utf-8",
    )

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