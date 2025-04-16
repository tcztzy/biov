"""DataFrame for biological data."""

from __future__ import annotations

from typing import Callable

from pandas import DataFrame

from .io.gff import GFFMixin


class BioDataFrame(GFFMixin, DataFrame):
    """DataFrame for biological data.

    Attributes:
        _gff_columns: Column names for GFF format. Length MUST be equal to 9.
    """

    @property
    def _constructor(self) -> Callable[..., BioDataFrame]:
        return BioDataFrame
