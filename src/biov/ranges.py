"""Range operations for BioDataFrame using PyRanges as a backend.

This module provides genomic range operations for BioDataFrame by leveraging
the PyRanges library. It includes conversion functions between BioDataFrame
and PyRanges objects, as well as wrapper methods for common range operations.
"""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Literal,
    Optional,
    Union,
)

import pyranges as pr
from pandas import DataFrame

if TYPE_CHECKING:
    from .dataframe import BioDataFrame


def to_pyranges(
    df: Union[DataFrame, "BioDataFrame"],
    seqid_col: str = "seqid",
    start_col: str = "start",
    end_col: str = "end",
    strand_col: Optional[str] = "strand",
) -> pr.PyRanges:
    """Convert a DataFrame or BioDataFrame to a PyRanges object.

    Args:
        df: DataFrame or BioDataFrame to convert
        seqid_col: Column name for sequence ID (chromosome)
        start_col: Column name for start position (0-based)
        end_col: Column name for end position (1-based)
        strand_col: Column name for strand information, or None if not available

    Returns:
        PyRanges object

    Notes:
        BioV uses 0-based start and 1-based end coordinates (BED-like),
        which is compatible with PyRanges' coordinate system.
    """
    # Create a copy to avoid modifying the original
    df_copy = df.copy()

    # Rename columns to match PyRanges expected format
    columns_map = {
        seqid_col: "Chromosome",
        start_col: "Start",
        end_col: "End",
    }

    if strand_col and strand_col in df.columns:
        columns_map[strand_col] = "Strand"

    df_copy = df_copy.rename(columns=columns_map)

    # Convert to PyRanges
    return pr.PyRanges(df_copy)


def from_pyranges(
    gr: pr.PyRanges,
    seqid_col: str = "seqid",
    start_col: str = "start",
    end_col: str = "end",
    strand_col: Optional[str] = "strand",
) -> "BioDataFrame":
    """Convert a PyRanges object to a BioDataFrame.

    Args:
        gr: PyRanges object to convert
        seqid_col: Column name for sequence ID (chromosome) in the output
        start_col: Column name for start position in the output (0-based)
        end_col: Column name for end position in the output (1-based)
        strand_col: Column name for strand information in the output, or None to exclude

    Returns:
        BioDataFrame with the converted data

    Notes:
        BioV uses 0-based start and 1-based end coordinates (BED-like),
        which is compatible with PyRanges' coordinate system.
    """
    # Convert PyRanges to DataFrame
    rdf = gr.as_df()

    # Rename columns to match BioDataFrame expected format
    columns_map = {
        "Chromosome": seqid_col,
        "Start": start_col,
        "End": end_col,
    }

    if strand_col and "Strand" in rdf.columns:
        columns_map["Strand"] = strand_col

    rdf = rdf.rename(columns=columns_map)

    # Convert to BioDataFrame
    from .dataframe import BioDataFrame

    return BioDataFrame(rdf)


class RangeMixin:
    """Mixin class providing genomic range operations for BioDataFrame.

    This mixin adds PyRanges-powered genomic range operations to BioDataFrame.
    """

    def to_pyranges(
        self,
        seqid_col: str = "seqid",
        start_col: str = "start",
        end_col: str = "end",
        strand_col: Optional[str] = "strand",
    ) -> pr.PyRanges:
        """Convert this BioDataFrame to a PyRanges object.

        Args:
            seqid_col: Column name for sequence ID (chromosome)
            start_col: Column name for start position (0-based)
            end_col: Column name for end position (1-based)
            strand_col: Column name for strand information, or None if not available

        Returns:
            PyRanges object
        """
        return to_pyranges(
            self,  # type: ignore
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )

    def overlap(
        self,
        other: Union["BioDataFrame", pr.PyRanges],
        how: Literal["first", "last", "containment", "member"] = "first",
        seqid_col: str = "seqid",
        start_col: str = "start",
        end_col: str = "end",
        strand_col: Optional[str] = "strand",
    ) -> "BioDataFrame":
        """Find overlapping ranges between this BioDataFrame and another.

        Args:
            other: Another BioDataFrame or PyRanges object to find overlaps with
            how: Method to use for overlap:
                - "first": Report first overlap (default)
                - "last": Report last overlap
                - "containment": Report ranges in self that contain ranges in other
                - "member": Report ranges in self that are contained in ranges in other
            seqid_col: Column name for sequence ID (chromosome)
            start_col: Column name for start position (0-based)
            end_col: Column name for end position (1-based)
            strand_col: Column name for strand information, or None if not available

        Returns:
            BioDataFrame with overlapping ranges
        """
        gr_self = self.to_pyranges(
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )

        if isinstance(other, pr.PyRanges):
            gr_other = other
        else:
            gr_other = other.to_pyranges(
                seqid_col=seqid_col,
                start_col=start_col,
                end_col=end_col,
                strand_col=strand_col,
            )

        result = gr_self.overlap(gr_other, how=how)

        return from_pyranges(
            result,
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )

    def intersect(
        self,
        other: Union["BioDataFrame", pr.PyRanges],
        seqid_col: str = "seqid",
        start_col: str = "start",
        end_col: str = "end",
        strand_col: Optional[str] = "strand",
    ) -> "BioDataFrame":
        """Find the intersection of ranges between this BioDataFrame and another.

        Args:
            other: Another BioDataFrame or PyRanges object to intersect with
            seqid_col: Column name for sequence ID (chromosome)
            start_col: Column name for start position (0-based)
            end_col: Column name for end position (1-based)
            strand_col: Column name for strand information, or None if not available

        Returns:
            BioDataFrame with intersecting ranges
        """
        gr_self = self.to_pyranges(
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )

        if isinstance(other, pr.PyRanges):
            gr_other = other
        else:
            gr_other = other.to_pyranges(
                seqid_col=seqid_col,
                start_col=start_col,
                end_col=end_col,
                strand_col=strand_col,
            )

        result = gr_self.intersect(gr_other)

        return from_pyranges(
            result,
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )

    def subtract_ranges(
        self,
        other: Union["BioDataFrame", pr.PyRanges],
        seqid_col: str = "seqid",
        start_col: str = "start",
        end_col: str = "end",
        strand_col: Optional[str] = "strand",
    ) -> "BioDataFrame":
        """Subtract ranges in other from this BioDataFrame.

        Args:
            other: Another BioDataFrame or PyRanges object to subtract
            seqid_col: Column name for sequence ID (chromosome)
            start_col: Column name for start position (0-based)
            end_col: Column name for end position (1-based)
            strand_col: Column name for strand information, or None if not available

        Returns:
            BioDataFrame with ranges after subtraction
        """
        gr_self = self.to_pyranges(
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )

        if isinstance(other, pr.PyRanges):
            gr_other = other
        else:
            gr_other = other.to_pyranges(
                seqid_col=seqid_col,
                start_col=start_col,
                end_col=end_col,
                strand_col=strand_col,
            )

        result = gr_self.subtract(gr_other)

        return from_pyranges(
            result,
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )

    def nearest(
        self,
        other: Union["BioDataFrame", pr.PyRanges],
        seqid_col: str = "seqid",
        start_col: str = "start",
        end_col: str = "end",
        strand_col: Optional[str] = "strand",
        suffix: str = "_b",
        how: Literal["upstream", "downstream", "next", "previous"] = "next",
    ) -> "BioDataFrame":
        """Find nearest ranges in other relative to this BioDataFrame.

        Args:
            other: Another BioDataFrame or PyRanges object to find nearest ranges in
            seqid_col: Column name for sequence ID (chromosome)
            start_col: Column name for start position (0-based)
            end_col: Column name for end position (1-based)
            strand_col: Column name for strand information, or None if not available
            suffix: Suffix to add to column names from the other DataFrame
            how: Method to use for finding nearest:
                - "upstream": Find nearest upstream range
                - "downstream": Find nearest downstream range
                - "next": Find nearest range in 3' direction (strand-aware)
                - "previous": Find nearest range in 5' direction (strand-aware)

        Returns:
            BioDataFrame with nearest ranges

        """
        gr_self = self.to_pyranges(
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )

        if isinstance(other, pr.PyRanges):
            gr_other = other
        else:
            gr_other = other.to_pyranges(
                seqid_col=seqid_col,
                start_col=start_col,
                end_col=end_col,
                strand_col=strand_col,
            )

        result = gr_self.nearest(gr_other, suffix=suffix, how=how)

        return from_pyranges(
            result,
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )
