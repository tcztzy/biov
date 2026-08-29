"""BioV-owned genomic range semantics backed by RuRanges NumPy kernels."""

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal, cast

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from pandas import DataFrame
from pandas.api.types import is_bool_dtype, is_integer_dtype
from ruranges import numpy as ruranges_numpy

if TYPE_CHECKING:
    from .dataframe import BioDataFrame

OverlapHow = Literal["first", "last", "containment", "member"]
NearestHow = Literal["upstream", "downstream", "next", "previous"]


@dataclass(frozen=True)
class _RangeArrays:
    """Validated arrays passed to RuRanges."""

    starts: NDArray[np.int64]
    ends: NDArray[np.int64]
    groups: NDArray[np.uint32]
    strands: NDArray[np.object_] | None


def _as_biodataframe(frame: DataFrame) -> "BioDataFrame":
    """Return a BioDataFrame with a fresh RangeIndex."""
    from .dataframe import BioDataFrame

    return BioDataFrame(frame.reset_index(drop=True))


def _validate_frame(
    frame: DataFrame,
    *,
    seqid_col: str,
    start_col: str,
    end_col: str,
) -> tuple[NDArray[np.int64], NDArray[np.int64], list[str]]:
    """Validate one half-open interval table and return kernel-ready values.

    Returns:
        Coordinate arrays and sequence IDs.

    Raises:
        ValueError: If required columns or valid interval values are missing.
    """
    missing = [
        column
        for column in (seqid_col, start_col, end_col)
        if column not in frame.columns
    ]
    if missing:
        raise ValueError(f"Missing range columns: {', '.join(missing)}")

    if frame.empty:
        return (
            np.empty(0, dtype=np.int64),
            np.empty(0, dtype=np.int64),
            [],
        )

    seqids = frame[seqid_col].tolist()
    if any(pd.isna(value) for value in seqids):
        raise ValueError(f"Range column {seqid_col!r} cannot contain missing values")
    if not all(isinstance(value, str) for value in seqids):
        raise ValueError(f"Range column {seqid_col!r} must contain strings")

    for column in (start_col, end_col):
        series = frame[column]
        if (
            not is_integer_dtype(series.dtype)
            or is_bool_dtype(series.dtype)
            or bool(series.isna().to_numpy().any())
        ):
            raise ValueError(f"Range column {column!r} must contain non-null integers")
        values = series.to_numpy(copy=False)
        if int(values.max()) > np.iinfo(np.int64).max:
            raise ValueError(
                f"Range column {column!r} exceeds signed 64-bit coordinates"
            )

    starts = frame[start_col].to_numpy(dtype=np.int64, copy=True)
    ends = frame[end_col].to_numpy(dtype=np.int64, copy=True)
    if np.any(starts < 0):
        raise ValueError("Range coordinates require 0 <= start")
    if np.any(starts >= ends):
        raise ValueError("Range coordinates require start < end")
    return starts, ends, seqids


def _encode_groups(
    left_seqids: list[str],
    right_seqids: list[str],
    left_strands: list[str] | None,
    right_strands: list[str] | None,
) -> tuple[NDArray[np.uint32], NDArray[np.uint32]]:
    """Encode exact chromosome or chromosome/strand keys across both operands.

    Returns:
        Shared integer group codes for the left and right operands.
    """
    mapping: dict[str | tuple[str, str], int] = {}

    def encode(seqids: list[str], strands: list[str] | None) -> NDArray[np.uint32]:
        result = np.empty(len(seqids), dtype=np.uint32)
        for index, seqid in enumerate(seqids):
            key: str | tuple[str, str]
            key = seqid if strands is None else (seqid, strands[index])
            code = mapping.setdefault(key, len(mapping))
            if code > np.iinfo(np.uint32).max:
                raise ValueError("Too many distinct range groups for RuRanges")
            result[index] = code
        return result

    return encode(left_seqids, left_strands), encode(right_seqids, right_strands)


def _prepare_operands(
    left: DataFrame,
    other: object,
    *,
    seqid_col: str,
    start_col: str,
    end_col: str,
    strand_col: str | None,
) -> tuple["BioDataFrame", _RangeArrays, _RangeArrays]:
    """Validate public operands and create shared group codes.

    Returns:
        Validated other frame and kernel arrays for both operands.

    Raises:
        TypeError: If other is not a BioDataFrame.
        ValueError: If columns, coordinates, or strand values are invalid.
    """
    from .dataframe import BioDataFrame

    if not isinstance(other, BioDataFrame):
        raise TypeError("Range operations require another BioDataFrame")

    left_starts, left_ends, left_seqids = _validate_frame(
        left,
        seqid_col=seqid_col,
        start_col=start_col,
        end_col=end_col,
    )
    right_starts, right_ends, right_seqids = _validate_frame(
        other,
        seqid_col=seqid_col,
        start_col=start_col,
        end_col=end_col,
    )

    left_strands: list[str] | None = None
    right_strands: list[str] | None = None
    if strand_col is not None:
        left_has_strand = strand_col in left.columns
        right_has_strand = strand_col in other.columns
        if left_has_strand != right_has_strand:
            raise ValueError(
                f"Selected strand column {strand_col!r} must exist on both operands or neither"
            )
        if left_has_strand:
            left_values = cast(list[str], left[strand_col].tolist())
            right_values = cast(list[str], other[strand_col].tolist())
            invalid = [
                value
                for value in (*left_values, *right_values)
                if value not in {"+", "-"}
            ]
            if invalid:
                raise ValueError("Range strand values must be '+' or '-'")
            left_strands = left_values
            right_strands = right_values

    left_groups, right_groups = _encode_groups(
        left_seqids,
        right_seqids,
        left_strands,
        right_strands,
    )
    left_strand_array = (
        None if left_strands is None else np.asarray(left_strands, dtype=object)
    )
    right_strand_array = (
        None if right_strands is None else np.asarray(right_strands, dtype=object)
    )
    return (
        other,
        _RangeArrays(left_starts, left_ends, left_groups, left_strand_array),
        _RangeArrays(right_starts, right_ends, right_groups, right_strand_array),
    )


def _overlap_pairs(
    left: _RangeArrays, right: _RangeArrays
) -> tuple[NDArray[np.uint32], NDArray[np.uint32]]:
    """Return all overlap pairs, bypassing RuRanges for empty inputs."""
    if left.starts.size == 0 or right.starts.size == 0:
        empty = np.empty(0, dtype=np.uint32)
        return empty, empty.copy()
    return cast(
        tuple[NDArray[np.uint32], NDArray[np.uint32]],
        ruranges_numpy.overlaps(
            starts=left.starts,
            ends=left.ends,
            starts2=right.starts,
            ends2=right.ends,
            groups=cast(Any, left.groups),
            groups2=cast(Any, right.groups),
            multiple="all",
            sort_output=True,
        ),
    )


def _ordered_pairs(
    left_indices: NDArray[np.uint32], right_indices: NDArray[np.uint32]
) -> tuple[NDArray[np.uint32], NDArray[np.uint32]]:
    """Order pairs by caller row and then other input row.

    Returns:
        Reordered left and right index arrays.
    """
    if left_indices.size == 0:
        return left_indices, right_indices
    order = np.lexsort((right_indices, left_indices))
    return left_indices[order], right_indices[order]


def _nearest_kernel(
    left: _RangeArrays,
    right: _RangeArrays,
    direction: Literal["forward", "backward"],
) -> tuple[NDArray[np.uint32], NDArray[np.uint32], NDArray[np.int64]]:
    """Return every tie at the nearest distance in one physical direction."""
    if left.starts.size == 0 or right.starts.size == 0:
        empty_index = np.empty(0, dtype=np.uint32)
        return empty_index, empty_index.copy(), np.empty(0, dtype=np.int64)
    return cast(
        tuple[NDArray[np.uint32], NDArray[np.uint32], NDArray[np.int64]],
        ruranges_numpy.nearest(
            starts=left.starts,
            ends=left.ends,
            starts2=right.starts,
            ends2=right.ends,
            groups=cast(Any, left.groups),
            groups2=cast(Any, right.groups),
            k=1,
            include_overlaps=True,
            direction=direction,
            ties="all",
            sort_output=True,
        ),
    )


def _choose_nearest(
    left_indices: NDArray[np.uint32],
    right_indices: NDArray[np.uint32],
    distances: NDArray[np.int64],
) -> tuple[NDArray[np.uint32], NDArray[np.uint32], NDArray[np.int64]]:
    """Choose the lowest other input row for each nearest-distance tie.

    Returns:
        One left index, right index, and distance per matched caller row.
    """
    if left_indices.size == 0:
        return left_indices, right_indices, distances
    order = np.lexsort((right_indices, distances, left_indices))
    left_indices = left_indices[order]
    right_indices = right_indices[order]
    distances = distances[order]
    keep = np.concatenate((np.array([True]), left_indices[1:] != left_indices[:-1]))
    return left_indices[keep], right_indices[keep], distances[keep]


def _nearest_result(
    left: DataFrame,
    right: DataFrame,
    left_indices: NDArray[np.uint32],
    right_indices: NDArray[np.uint32],
    distances: NDArray[np.int64],
    *,
    seqid_col: str,
    suffix: str,
) -> "BioDataFrame":
    """Append selected other columns with a deterministic collision policy.

    Returns:
        Joined nearest-neighbor BioDataFrame.

    Raises:
        ValueError: If suffixing cannot produce unique output columns.
    """
    result = left.iloc[left_indices].reset_index(drop=True).copy()
    selected_right = right.iloc[right_indices].reset_index(drop=True)
    targets: list[tuple[str, str]] = []
    occupied = set(result.columns)
    for column in right.columns:
        if column == seqid_col:
            continue
        target = f"{column}{suffix}" if column in left.columns else str(column)
        if target == "Distance" or target in occupied:
            raise ValueError(
                f"Nearest output column {target!r} is not unique; choose another suffix"
            )
        occupied.add(target)
        targets.append((str(column), target))

    for source, target in targets:
        result[target] = selected_right[source]
    if distances.size:
        result["Distance"] = distances.astype(np.int64, copy=False)
    else:
        result["Distance"] = pd.Series(dtype="int64")
    return _as_biodataframe(result)


class RangeMixin:
    """Genomic range operations for :class:`~biov.dataframe.BioDataFrame`."""

    def overlap(
        self,
        other: "BioDataFrame",
        how: OverlapHow = "first",
        seqid_col: str = "seqid",
        start_col: str = "start",
        end_col: str = "end",
        strand_col: str | None = "strand",
    ) -> "BioDataFrame":
        """Return each caller row matching the requested overlap mode once.

        Raises:
            ValueError: If the overlap mode or range operands are invalid.
        """
        if how not in {"first", "last", "containment", "member"}:
            raise ValueError(f"Unsupported overlap mode: {how!r}")
        frame = cast(DataFrame, self)
        other, left, right = _prepare_operands(
            frame,
            other,
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )
        left_indices, right_indices = _overlap_pairs(left, right)
        if how == "containment":
            keep = (left.starts[left_indices] <= right.starts[right_indices]) & (
                left.ends[left_indices] >= right.ends[right_indices]
            )
            left_indices = left_indices[keep]
        elif how == "member":
            keep = (left.starts[left_indices] >= right.starts[right_indices]) & (
                left.ends[left_indices] <= right.ends[right_indices]
            )
            left_indices = left_indices[keep]
        selected = np.unique(left_indices)
        return _as_biodataframe(frame.iloc[selected].copy())

    def intersect(
        self,
        other: "BioDataFrame",
        seqid_col: str = "seqid",
        start_col: str = "start",
        end_col: str = "end",
        strand_col: str | None = "strand",
    ) -> "BioDataFrame":
        """Return one clipped caller-metadata row per overlap pair."""
        frame = cast(DataFrame, self)
        other, left, right = _prepare_operands(
            frame,
            other,
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )
        left_indices, right_indices = _ordered_pairs(*_overlap_pairs(left, right))
        if left_indices.size == 0:
            return _as_biodataframe(frame.iloc[:0].copy())
        result = frame.iloc[left_indices].reset_index(drop=True).copy()
        result[start_col] = np.maximum(
            left.starts[left_indices], right.starts[right_indices]
        )
        result[end_col] = np.minimum(left.ends[left_indices], right.ends[right_indices])
        return _as_biodataframe(result)

    def subtract_ranges(
        self,
        other: "BioDataFrame",
        seqid_col: str = "seqid",
        start_col: str = "start",
        end_col: str = "end",
        strand_col: str | None = "strand",
    ) -> "BioDataFrame":
        """Subtract matching other intervals and return all residual fragments.

        Returns:
            Residual interval fragments with caller metadata.
        """
        frame = cast(DataFrame, self)
        other, left, right = _prepare_operands(
            frame,
            other,
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )
        if left.starts.size == 0:
            return _as_biodataframe(frame.iloc[:0].copy())
        if right.starts.size == 0:
            return _as_biodataframe(frame.copy())
        indices, starts, ends = cast(
            tuple[NDArray[np.uint32], NDArray[np.int64], NDArray[np.int64]],
            ruranges_numpy.subtract(
                left.starts,
                left.ends,
                right.starts,
                right.ends,
                groups=cast(Any, left.groups),
                groups2=cast(Any, right.groups),
                sort_output=True,
            ),
        )
        order = np.lexsort((starts, indices))
        indices, starts, ends = indices[order], starts[order], ends[order]
        result = frame.iloc[indices].reset_index(drop=True).copy()
        result[start_col] = starts
        result[end_col] = ends
        return _as_biodataframe(result)

    def nearest(
        self,
        other: "BioDataFrame",
        seqid_col: str = "seqid",
        start_col: str = "start",
        end_col: str = "end",
        strand_col: str | None = "strand",
        suffix: str = "_b",
        how: NearestHow = "next",
    ) -> "BioDataFrame":
        """Append one nearest other row under explicit physical/strand directions.

        Returns:
            Caller rows joined to at most one nearest other row.

        Raises:
            ValueError: If direction, strands, range values, or suffix are invalid.
        """
        if how not in {"upstream", "downstream", "next", "previous"}:
            raise ValueError(f"Unsupported nearest direction: {how!r}")
        frame = cast(DataFrame, self)
        other, left, right = _prepare_operands(
            frame,
            other,
            seqid_col=seqid_col,
            start_col=start_col,
            end_col=end_col,
            strand_col=strand_col,
        )

        if how in {"next", "previous"}:
            direction: Literal["forward", "backward"]
            direction = "forward" if how == "next" else "backward"
            left_indices, right_indices, distances = _nearest_kernel(
                left, right, direction
            )
        else:
            if left.strands is None or right.strands is None:
                raise ValueError(
                    "Nearest upstream/downstream directions require strand columns"
                )
            forward = _nearest_kernel(left, right, "forward")
            backward = _nearest_kernel(left, right, "backward")
            plus_uses_forward = how == "downstream"
            forward_keep = (left.strands[forward[0]] == "+") == plus_uses_forward
            backward_keep = (left.strands[backward[0]] == "+") != plus_uses_forward
            left_indices = np.concatenate(
                (forward[0][forward_keep], backward[0][backward_keep])
            )
            right_indices = np.concatenate(
                (forward[1][forward_keep], backward[1][backward_keep])
            )
            distances = np.concatenate(
                (forward[2][forward_keep], backward[2][backward_keep])
            )

        left_indices, right_indices, distances = _choose_nearest(
            left_indices, right_indices, distances
        )
        return _nearest_result(
            frame,
            other,
            left_indices,
            right_indices,
            distances,
            seqid_col=seqid_col,
            suffix=suffix,
        )
