"""Refresh the checkout's registry asset through the BioV CLI."""

import sys
from pathlib import Path

from biov.cli import app
from biov.registry import REGISTRY_ASSET_NAME

DEFAULT_OUTPUT = (
    Path(__file__).parents[1] / "src" / "biov" / "assets" / REGISTRY_ASSET_NAME
)


def main() -> None:
    """Run the registry asset updater."""
    app(["update-identifiers-registry", "--output", str(DEFAULT_OUTPUT), *sys.argv[1:]])


if __name__ == "__main__":
    main()
