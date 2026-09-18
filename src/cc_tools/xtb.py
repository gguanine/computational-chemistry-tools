from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from collections.abc import Sequence

from .molecule import Molecule
from .xyz import read_xyz, write_xyz


def optimize(
    molecule: Molecule,
    *,
    method: str = "gfn2",
    extra_args: Sequence[str] = (),
    executable: str = "xtb",
) -> Molecule:
    path_xtb = shutil.which(executable)

    if path_xtb is None:
        raise FileNotFoundError(
            f"xTB executable {executable!r} not found in PATH"
        )

    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        input_path = workdir / "input.xyz"
        output_path = workdir / "xtbopt.xyz"

        write_xyz(molecule, input_path)

        command = [
            path_xtb,
            input_path.name,
            "--opt",
            "--gfn",
            method.removeprefix("gfn"),
        ]

        if molecule.charge != 0:
            command += ["--chrg", str(molecule.charge)]

        uhf = molecule.multiplicity - 1
        if uhf != 0:
            command += ["--uhf", str(uhf)]

        command += list(extra_args)

        try:
            completed = subprocess.run(
                command,
                cwd=workdir,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                "xTB optimization failed\n"
                f"command: {' '.join(command)}\n"
                f"stdout:\n{exc.stdout}\n"
                f"stderr:\n{exc.stderr}"
            ) from exc

        if not output_path.exists():
            raise RuntimeError(
                "xTB finished successfully but xtbopt.xyz was not produced.\n"
                f"stdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}"
            )

        optimized = read_xyz(
            output_path,
            charge=molecule.charge,
            multiplicity=molecule.multiplicity,
            name=molecule.name,
            smiles=molecule.smiles,
        )

        return optimized