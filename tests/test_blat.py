"""Tests for the BLAT command-line interface."""

from typing import Any, cast

from typer.main import get_command

from biov.executables import blat_app


def test_blat_literal_options_expose_exact_choices() -> None:
    """Derive all finite BLAT option domains from Literal annotations."""
    expected = {
        "t": ("dna", "prot", "dnax"),
        "q": ("dna", "rna", "prot", "dnax", "rnax"),
        "oneOff": ("0", "1"),
        "mask": ("lower", "upper", "out", "file.out"),
        "qMask": ("lower", "upper", "out", "file.out"),
        "repeats": ("lower", "upper", "out", "file.out"),
        "out": (
            "psl",
            "pslx",
            "axt",
            "maf",
            "sim4",
            "wublast",
            "blast",
            "blast8",
            "blast9",
        ),
    }
    parameters = {
        parameter.name: parameter for parameter in get_command(blat_app).params
    }

    for name, choices in expected.items():
        parameter_type = cast(Any, parameters[name].type)
        assert tuple(parameter_type.choices) == choices
