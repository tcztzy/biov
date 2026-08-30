"""DataFrame for biological data."""

from pandas import DataFrame

from .io.gff import GFFMixin
from .ranges import RangeMixin


# pandas 3 narrows NDFrame's Self returns to DataFrame in these overrides.
class BioDataFrame(GFFMixin, RangeMixin, DataFrame):  # ty: ignore[invalid-method-override]
    """DataFrame for biological data.

    Attributes:
        _gff_columns: Column names for GFF format. Length MUST be equal to 9.
    """

    _metadata = ["_gff_columns"]  # noqa: RUF012

    @property
    def _constructor(self) -> type["BioDataFrame"]:
        return BioDataFrame
