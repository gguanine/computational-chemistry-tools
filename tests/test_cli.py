from pathlib import Path

from cc_tools.cli import _build_parser


def test_gaussian_prepare_sp_parser():
    parser = _build_parser()

    args = parser.parse_args([
        "gaussian",
        "prepare-sp",
        "opt.log",
        "sp.gjf",
        "--route",
        "#p B3LYP/6-31G(d)",
        "--require-minimum",
        "--nprocs",
        "8",
        "--memory",
        "16GB",
    ])

    assert args.command == "gaussian"
    assert args.gaussian_command == "prepare-sp"

    assert args.input == Path("opt.log")
    assert args.output == Path("sp.gjf")

    assert args.route == "#p B3LYP/6-31G(d)"
    assert args.require_minimum is True
    assert args.nprocs == 8
    assert args.memory == "16GB"

