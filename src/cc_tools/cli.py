from __future__ import annotations

import argparse
from pathlib import Path

from .gaussian import prepare_single_point, prepare_single_points, read_log


def _cmd_gaussian_prepare_sp(args: argparse.Namespace) -> None:
    if args.input.is_dir():
        prepare_single_points(
            args.input,
            args.output,
            route=args.route,
            require_minimum=args.require_minimum,
            nprocs=args.nprocs,
            memory=args.memory,
            checkpoint=args.checkpoint,
        )
        return

    result = read_log(args.input)

    checkpoint_name = (
        args.output.with_suffix(".chk").name
        if args.checkpoint
        else None
    )

    gaussian_input = prepare_single_point(
        result,
        route=args.route,
        require_minimum=args.require_minimum,
        nprocs=args.nprocs,
        memory=args.memory,
        checkpoint=checkpoint_name,
    )

    gaussian_input.write(args.output)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cc-tools",
        description="Lightweight computational chemistry tools.",
    )

    commands = parser.add_subparsers(
        dest="command",
        required=True,
    )

    # ---------------------------------------------------------
    # gaussian
    # ---------------------------------------------------------

    gaussian_parser = commands.add_parser(
        "gaussian",
        help="Gaussian utilities.",
    )

    gaussian_commands = gaussian_parser.add_subparsers(
        dest="gaussian_command",
        required=True,
    )

    # ---------------------------------------------------------
    # gaussian prepare-sp
    # ---------------------------------------------------------

    prepare_sp_parser = gaussian_commands.add_parser(
        "prepare-sp",
        help="Prepare a Gaussian single-point input from an Opt/Freq log.",
    )

    prepare_sp_parser.add_argument(
        "input",
        type=Path,
        help="Gaussian Opt/Freq log file.",
    )

    prepare_sp_parser.add_argument(
        "output",
        type=Path,
        help="Output Gaussian input file.",
    )

    prepare_sp_parser.add_argument(
        "--route",
        required=True,
        help="Gaussian route section.",
    )

    prepare_sp_parser.add_argument(
        "--require-minimum",
        action="store_true",
        help="Require no imaginary frequencies.",
    )

    prepare_sp_parser.add_argument(
        "--nprocs",
        type=int,
        help="Number of shared-memory processors.",
    )

    prepare_sp_parser.add_argument(
        "--memory",
        help="Gaussian memory specification, e.g. 8GB.",
    )

    prepare_sp_parser.add_argument(
        "--checkpoint",
        help="Gaussian checkpoint file name.",
    )

    prepare_sp_parser.set_defaults(
        func=_cmd_gaussian_prepare_sp,
    )

    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    try:
        args.func(args)
    except (RuntimeError, ValueError, FileNotFoundError) as exc:
        parser.exit(
            status=1,
            message=f"error: {exc}\n",
        )


if __name__ == "__main__":
    main()