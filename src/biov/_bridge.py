"""Paired-interpreter compatibility boundary for the native BioV executable.

This private module is invoked only by the installed Rust command. It preserves
legacy Typer argv and stdio while rejecting an interpreter that imports a different
BioV distribution. Python capabilities remain implemented in ``biov.cli``.
"""

import sys
from importlib.metadata import PackageNotFoundError, distribution
from pathlib import Path


def validate_pairing(binary: Path) -> None:
    """Require the calling executable to belong to this installed distribution.

    Args:
        binary: Real native executable path supplied by the Rust router.

    Raises:
        ValueError: BioV metadata is missing or records a different executable.
    """
    try:
        installed = distribution("biov")
    except PackageNotFoundError as error:
        raise ValueError(
            "paired Python interpreter has no installed BioV distribution"
        ) from error
    expected = binary.resolve(strict=True)
    for entry in installed.files or ():
        if entry.name not in {"biov", "biov.exe"}:
            continue
        candidate = installed.locate_file(entry)
        try:
            if candidate.resolve(strict=True) == expected:
                return
        except OSError:
            continue
    raise ValueError(
        "paired Python interpreter does not own this BioV executable; "
        "reinstall BioV with `uv tool install --force biov`"
    )


def main() -> None:
    """Verify pairing and dispatch unchanged legacy arguments to Typer.

    Raises:
        SystemExit: With status 2 when invocation or distribution pairing fails.
    """
    if not sys.argv or Path(sys.argv[0]).name not in {"biov", "biov.exe"}:
        sys.stderr.write(
            "biov: Python bridge requires the installed native executable\n"
        )
        raise SystemExit(2)
    try:
        validate_pairing(Path(sys.argv[0]))
    except (OSError, ValueError) as error:
        sys.stderr.write(f"biov: {error}\n")
        raise SystemExit(2) from error
    from .cli import app

    app(args=sys.argv[1:], prog_name="biov")


if __name__ == "__main__":
    main()
