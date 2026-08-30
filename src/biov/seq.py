"""Scalar and explicitly typed pandas sequence APIs for BioV."""

import builtins
from collections.abc import Sequence
from typing import Any, Literal, Self, cast, overload

import numpy as np
import numpy.typing as npt
import pandas as pd
from Bio.Seq import Seq as _Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from numpy.typing import NDArray
from pandas import DataFrame, Series
from pandas.api.extensions import (
    ExtensionArray,
    ExtensionDtype,
    register_extension_dtype,
    register_series_accessor,
    take,
)
from pandas.api.types import is_scalar, pandas_dtype
from pandas.api.typing.aliases import (
    ArrayLike,
    AstypeArg,
    Dtype,
    ScalarIndexer,
    SequenceIndexer,
    TakeIndexer,
)
from pydantic import GetCoreSchemaHandler
from pydantic_core import CoreSchema, core_schema

SequenceKind = Literal["dna", "rna", "protein"]

_DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
_RNA_ALPHABET = frozenset("ACGURYSWKMBDHVN")
_CANONICAL_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
_PROTEIN_ALPHABET = frozenset(f"{_CANONICAL_AMINO_ACIDS}BJOUXZ*")
_ALPHABETS = {
    "dna": _DNA_ALPHABET,
    "rna": _RNA_ALPHABET,
    "protein": _PROTEIN_ALPHABET,
}


class SequenceValidationError(ValueError):
    """Raised when a value violates a declared BioV sequence contract."""


class Seq(_Seq):
    """Biopython sequence with Pydantic validation and serialization."""

    @classmethod
    def __get_pydantic_core_schema__(
        cls, source_type: Any, handler: GetCoreSchemaHandler
    ) -> CoreSchema:
        """Return the Pydantic validation and string serialization schema."""
        return core_schema.no_info_plain_validator_function(
            _Seq,
            serialization=core_schema.to_string_ser_schema(),
        )


@register_extension_dtype
class SequenceDtype(ExtensionDtype):
    """Pandas dtype carrying one explicit DNA, RNA, or protein kind."""

    _metadata = ("sequence_kind",)

    def __init__(self, sequence_kind: SequenceKind) -> None:
        """Create one explicit DNA, RNA, or protein dtype.

        Raises:
            TypeError: If sequence_kind is unknown.
        """
        if sequence_kind not in _ALPHABETS:
            raise TypeError(f"Unknown BioV sequence kind: {sequence_kind!r}")
        self.sequence_kind = sequence_kind

    @property
    def name(self) -> str:
        """Registered pandas dtype name."""
        return f"biov.{self.sequence_kind}"

    @property
    def type(self) -> builtins.type[str]:
        """Scalar Python type."""
        return str

    @property
    def kind(self) -> str:
        """NumPy object-kind marker."""
        return "O"

    @property
    def na_value(self) -> pd.api.typing.NAType:
        """Pandas scalar missing value."""
        return pd.NA

    @classmethod
    def construct_from_string(cls, string: str) -> Self:
        """Construct a dtype from its exact registered string name.

        Returns:
            Matching BioV sequence dtype.

        Raises:
            TypeError: If string is not a registered BioV sequence dtype.
        """
        prefix = "biov."
        if not isinstance(string, str) or not string.startswith(prefix):
            raise TypeError(f"Cannot construct a SequenceDtype from {string!r}")
        kind = string.removeprefix(prefix)
        if kind not in _ALPHABETS:
            raise TypeError(f"Cannot construct a SequenceDtype from {string!r}")
        return cls(cast(SequenceKind, kind))

    def construct_array_type(self) -> builtins.type["SequenceArray"]:
        """Return the ExtensionArray paired with this dtype."""
        return SequenceArray


def _coerce_dtype(dtype: object) -> SequenceDtype:
    """Return one BioV dtype or reject an untyped array construction.

    Raises:
        TypeError: If dtype does not declare a BioV sequence kind.
    """
    if isinstance(dtype, SequenceDtype):
        return dtype
    if isinstance(dtype, str):
        return SequenceDtype.construct_from_string(dtype)
    raise TypeError("BioV sequence arrays require an explicit sequence dtype")


def _is_missing(value: object) -> bool:
    """Return whether one scalar is a pandas-compatible missing value."""
    if value is None or value is pd.NA:
        return True
    missing = pd.isna(value)
    return isinstance(missing, (bool, np.bool_)) and bool(missing)


def _normalize_sequence(
    value: object, dtype: SequenceDtype
) -> str | pd.api.typing.NAType:
    """Normalize one scalar under its declared alphabet.

    Returns:
        Uppercase sequence or pandas missing value.

    Raises:
        SequenceValidationError: If value is not a string or violates its alphabet.
    """
    if _is_missing(value):
        return pd.NA
    if not isinstance(value, str):
        raise SequenceValidationError(
            f"{dtype.name} values must be a string or missing, got {type(value).__name__}"
        )
    normalized = value.upper()
    invalid = sorted(set(normalized) - _ALPHABETS[dtype.sequence_kind])
    if invalid:
        label = (
            dtype.sequence_kind.upper()
            if dtype.sequence_kind in {"dna", "rna"}
            else "protein"
        )
        raise SequenceValidationError(
            f"invalid {label} symbol(s): {', '.join(invalid)}"
        )
    return normalized


class SequenceArray(ExtensionArray):
    """Minimal nullable sequence array storing normalized Python strings."""

    def __init__(
        self,
        values: Sequence[object] | NDArray[np.object_],
        dtype: SequenceDtype,
        *,
        copy: bool = False,
    ) -> None:
        """Validate and store normalized sequence values."""
        self._dtype = _coerce_dtype(dtype)
        normalized = [_normalize_sequence(value, self._dtype) for value in values]
        self._data = np.asarray(normalized, dtype=object)
        if copy:
            self._data = self._data.copy()

    @classmethod
    def _from_sequence(
        cls,
        scalars: Sequence[object],
        *,
        dtype: Dtype | None = None,
        copy: bool = False,
    ) -> Self:
        """Construct from pandas scalars under an explicit BioV dtype.

        Returns:
            Validated sequence array.
        """
        return cls(scalars, _coerce_dtype(dtype), copy=copy)

    @classmethod
    def _from_sequence_of_strings(
        cls,
        strings: Sequence[str],
        *,
        dtype: ExtensionDtype,
        copy: bool = False,
    ) -> Self:
        """Construct from parser strings under an explicit BioV dtype.

        Returns:
            Validated sequence array.
        """
        return cls(strings, _coerce_dtype(dtype), copy=copy)

    @classmethod
    def _from_factorized(
        cls, values: NDArray[np.object_], original: ExtensionArray
    ) -> Self:
        """Reconstruct an array after pandas factorization.

        Returns:
            Reconstructed sequence array.

        Raises:
            TypeError: If original is not a SequenceArray.
        """
        if not isinstance(original, SequenceArray):
            raise TypeError("Sequence factorization requires a SequenceArray")
        return cls(values, original.dtype)

    @classmethod
    def _concat_same_type(cls, to_concat: Sequence[ExtensionArray]) -> Self:
        """Concatenate arrays while preserving one exact sequence dtype.

        Returns:
            Concatenated sequence array.

        Raises:
            TypeError: If inputs are empty, untyped, or have differing dtypes.
        """
        arrays = [array for array in to_concat if isinstance(array, SequenceArray)]
        if len(arrays) != len(to_concat) or not arrays:
            raise TypeError("Sequence concatenation requires SequenceArray inputs")
        dtype = arrays[0].dtype
        if any(array.dtype != dtype for array in arrays[1:]):
            raise TypeError("Cannot concatenate different BioV sequence dtypes")
        return cls(np.concatenate([array._data for array in arrays]), dtype)

    @property
    def dtype(self) -> SequenceDtype:
        """Explicit sequence dtype for this array."""
        return self._dtype

    @property
    def _can_hold_na(self) -> bool:
        """Whether this array supports missing values."""
        return True

    @property
    def nbytes(self) -> int:
        """Bytes used by the object-reference array."""
        return self._data.nbytes

    def __len__(self) -> int:
        """Return array length."""
        return len(self._data)

    @overload
    def __getitem__(self, item: ScalarIndexer) -> Any: ...

    @overload
    def __getitem__(  # type: ignore[overload-cannot-match]
        self, item: SequenceIndexer
    ) -> Self: ...

    def __getitem__(self, item: ScalarIndexer | SequenceIndexer) -> Self | Any:
        """Return one scalar or a sliced sequence array."""
        result = self._data[item]
        if isinstance(item, (int, np.integer)):
            return result
        array = type(self)(result, self.dtype)
        array._readonly = self._readonly
        return array

    def __setitem__(self, key: object, value: object) -> None:
        """Validate and assign one or more sequence values.

        Raises:
            ValueError: If pandas exposed this array as read-only.
        """
        if self._readonly:
            raise ValueError("Cannot modify read-only array")
        if is_scalar(value):
            self._data[cast(Any, key)] = _normalize_sequence(value, self.dtype)
            return
        values = [
            _normalize_sequence(item, self.dtype)
            for item in cast(Sequence[object], value)
        ]
        self._data[cast(Any, key)] = values

    def __array__(
        self, dtype: np.dtype[Any] | None = None, copy: bool | None = None
    ) -> NDArray[Any]:
        """Return a NumPy object representation."""
        if copy is True:
            return np.array(self._data, dtype=dtype, copy=True)
        result = np.asarray(self._data, dtype=dtype)
        if self._readonly and np.shares_memory(result, self._data):
            result = result.view()
            result.flags.writeable = False
        return result

    def __eq__(self, other: object) -> Any:
        """Compare sequence values with nullable string semantics.

        Returns:
            Nullable elementwise equality values.
        """
        left = pd.array(self._data, dtype="string[python]")
        right = (
            pd.array(other._data, dtype="string[python]")
            if isinstance(other, SequenceArray)
            else other
        )
        return left == right

    def __ne__(self, other: object) -> Any:
        """Compare sequence inequality with nullable string semantics.

        Returns:
            Nullable elementwise inequality values.
        """
        left = pd.array(self._data, dtype="string[python]")
        right = (
            pd.array(other._data, dtype="string[python]")
            if isinstance(other, SequenceArray)
            else other
        )
        return left != right

    def isna(self) -> NDArray[np.bool_]:
        """Return the missing-value mask."""
        return np.fromiter(
            (value is pd.NA for value in self._data),
            dtype=np.bool_,
            count=len(self),
        )

    def take(
        self,
        indices: TakeIndexer,
        *,
        allow_fill: bool = False,
        fill_value: object = None,
    ) -> Self:
        """Take positional values using pandas' fill convention.

        Returns:
            Sequence array containing requested positions.
        """
        if allow_fill and fill_value is None:
            fill_value = self.dtype.na_value
        values = take(
            self._data,
            indices,
            allow_fill=allow_fill,
            fill_value=fill_value,
        )
        return type(self)(cast(NDArray[np.object_], values), self.dtype)

    def copy(self) -> Self:
        """Return an independent sequence array."""
        return type(self)(self._data.copy(), self.dtype)

    @overload
    def astype(self, dtype: npt.DTypeLike, copy: bool = True) -> NDArray[Any]: ...

    @overload
    def astype(  # type: ignore[overload-cannot-match]
        self, dtype: ExtensionDtype, copy: bool = True
    ) -> ExtensionArray: ...

    @overload
    def astype(  # type: ignore[overload-cannot-match]
        self, dtype: AstypeArg, copy: bool = True
    ) -> ArrayLike: ...

    def astype(self, dtype: AstypeArg, copy: bool = True) -> ArrayLike:
        """Cast to another sequence dtype or a NumPy-compatible dtype.

        Returns:
            SequenceArray or NumPy array with requested dtype.
        """
        resolved = pandas_dtype(dtype)
        if isinstance(resolved, SequenceDtype):
            if resolved == self.dtype:
                return self.copy() if copy else self
            return type(self)(self._data, resolved, copy=copy)
        if isinstance(resolved, ExtensionDtype):
            return pd.array(self._data, dtype=resolved, copy=copy)
        return np.asarray(self._data).astype(cast(Any, resolved), copy=copy)

    def _values_for_factorize(self) -> tuple[NDArray[np.object_], object]:
        values = self._data.copy()
        values[self.isna()] = None
        return values, None


@register_series_accessor("seq")
class SequenceAccessor:
    """Vectorized Biopython algorithms for explicitly typed sequence Series."""

    def __init__(self, pandas_obj: Series) -> None:
        """Validate and retain an explicitly typed sequence Series.

        Raises:
            AttributeError: If Series lacks an explicit BioV sequence dtype.
        """
        if not isinstance(pandas_obj.dtype, SequenceDtype) or not isinstance(
            pandas_obj.array, SequenceArray
        ):
            raise AttributeError(  # noqa: TRY004 -- pandas accessor protocol.
                ".seq requires a BioV sequence dtype: biov.dna, biov.rna, or biov.protein"
            )
        self._obj = pandas_obj
        self._array = pandas_obj.array

    @property
    def length(self) -> Series:
        """Nullable sequence lengths."""
        values = [
            pd.NA if value is pd.NA else len(cast(str, value)) for value in self._array
        ]
        return Series(values, index=self._obj.index, dtype="Int64", name=self._obj.name)

    def _require_nucleic(self) -> None:
        if self._array.dtype.sequence_kind not in {"dna", "rna"}:
            raise TypeError("This .seq operation requires a DNA or RNA dtype")

    def _require_protein(self) -> None:
        if self._array.dtype.sequence_kind != "protein":
            raise TypeError("This .seq operation requires a protein dtype")

    def _validated_proteins(self) -> list[str | None]:
        """Return canonical proteins or one missing marker per row.

        Raises:
            SequenceValidationError: If a present protein is empty or noncanonical.
        """
        self._require_protein()
        result: list[str | None] = []
        canonical = frozenset(_CANONICAL_AMINO_ACIDS)
        for value in self._array:
            if value is pd.NA:
                result.append(None)
                continue
            sequence = cast(str, value)
            if not sequence or not set(sequence) <= canonical:
                raise SequenceValidationError(
                    "Protein analysis requires a non-empty canonical amino-acid sequence"
                )
            result.append(sequence)
        return result

    def reverse_complement(self) -> Series:
        """Return DNA or RNA reverse complements with the same dtype."""
        self._require_nucleic()
        is_rna = self._array.dtype.sequence_kind == "rna"
        values = [
            pd.NA
            if value is pd.NA
            else str(
                _Seq(cast(str, value)).reverse_complement_rna()
                if is_rna
                else _Seq(cast(str, value)).reverse_complement()
            )
            for value in self._array
        ]
        return Series(
            values,
            index=self._obj.index,
            dtype=self._array.dtype,
            name=self._obj.name,
        )

    def gc_fraction(self) -> Series:
        """Return Biopython weighted-IUPAC GC fractions."""
        self._require_nucleic()
        values = [
            pd.NA
            if value is pd.NA
            else gc_fraction(cast(str, value), ambiguous="weighted")
            for value in self._array
        ]
        return Series(
            values, index=self._obj.index, dtype="Float64", name=self._obj.name
        )

    def translate(self, table: int | str = 1, to_stop: bool = False) -> Series:
        """Translate complete DNA/RNA codons into typed protein sequences.

        Returns:
            Nullable protein sequence Series.

        Raises:
            SequenceValidationError: If a present sequence has an incomplete codon.
        """
        self._require_nucleic()
        values: list[str | pd.api.typing.NAType] = []
        for value in self._array:
            if value is pd.NA:
                values.append(pd.NA)
                continue
            sequence = cast(str, value)
            if len(sequence) % 3:
                raise SequenceValidationError(
                    "Nucleic-acid translation requires complete codons"
                )
            values.append(
                str(_Seq(sequence).translate(table=cast(Any, table), to_stop=to_stop))
            )
        return Series(
            values,
            index=self._obj.index,
            dtype=SequenceDtype("protein"),
            name=self._obj.name,
        )

    def molecular_weight(self) -> Series:
        """Return Biopython molecular weights for canonical proteins."""
        values = [
            pd.NA if sequence is None else ProteinAnalysis(sequence).molecular_weight()
            for sequence in self._validated_proteins()
        ]
        return Series(
            values, index=self._obj.index, dtype="Float64", name=self._obj.name
        )

    def isoelectric_point(self) -> Series:
        """Return Biopython isoelectric points for canonical proteins."""
        values = [
            pd.NA if sequence is None else ProteinAnalysis(sequence).isoelectric_point()
            for sequence in self._validated_proteins()
        ]
        return Series(
            values, index=self._obj.index, dtype="Float64", name=self._obj.name
        )

    def amino_acid_composition(self) -> DataFrame:
        """Return canonical amino-acid percentages in stable column order."""
        rows: list[dict[str, float | pd.api.typing.NAType]] = []
        for sequence in self._validated_proteins():
            if sequence is None:
                rows.append(dict.fromkeys(_CANONICAL_AMINO_ACIDS, pd.NA))
                continue
            percentages = ProteinAnalysis(sequence).amino_acids_percent
            rows.append(
                {
                    amino_acid: float(percentages[amino_acid])
                    for amino_acid in _CANONICAL_AMINO_ACIDS
                }
            )
        return DataFrame(
            rows,
            columns=pd.Index(list(_CANONICAL_AMINO_ACIDS)),
            index=self._obj.index,
            dtype="Float64",
        )
