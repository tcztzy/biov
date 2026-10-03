"""Command-line entry point for the BioV-owned crisprprimer implementation."""

from pathlib import Path
from typing import Annotated

import typer

from . import crisprprimer
from .nuclease import Nucleases

app = typer.Typer()


# MOTIVATION: why defined main function instead of directly decorating the crisprprimer function?
# ANSWER: need use crisprprimer function in other places.
@app.command()
def main(
    output: Path,
    regions: Annotated[
        list[str] | None,
        typer.Option(
            "-r",
            help="The region to search for CRISPR primers, Locus_id (in GFF's attributes) is also accepted. Format: 'chr02:21250083..21251343' or 'Os02g0559800' or 'LOC_Os02g35329'.",
        ),
    ] = None,
    system: Annotated[
        Nucleases,
        typer.Option(
            help="The system of the sequence.",
            show_choices=True,
        ),
    ] = Nucleases.SpCas9,
) -> None:
    """Run the command-line workflow with the caller-selected inputs.

    Raises:
        ValueError: If the output has an unsupported suffix.
    """
    match output.suffix:
        case "":
            output.mkdir(parents=True, exist_ok=True)
            crisprprimer(regions=regions or None, system=system, blocks_dir=output)
        case _:
            raise ValueError(f"Suffix {output.suffix} is not supported!")


def run() -> None:
    """Invoke the crisprprimer command-line application."""
    app()


if __name__ == "__main__":
    run()
