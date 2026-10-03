"""Independent scientific and boundary checks for native sequence metrics.

IUPAC meanings: INSDC Feature Table 11.4 (April 2026), section 7.4.1,
https://www.insdc.org/submitting-standards/feature-table/ . GC expectations are
rational probabilities derived from base sets, not an implementation weight
lookup. Biopython 1.88 supplies the pinned differential reference in uv.lock.
See docs/guides/sequence-contract.md for semantics and distribution gates.
"""

import importlib
import inspect
import math
import random
from collections.abc import Callable, Sequence
from fractions import Fraction
from importlib.machinery import EXTENSION_SUFFIXES
from itertools import product
from pathlib import Path

import Bio
import Bio.SeqUtils
import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq as ReferenceSeq

import biov

REFERENCE_BIOPYTHON_VERSION = "1.88"
DNA_ALPHABET = "ACGTRYSWKMBDHVN"
RNA_ALPHABET = "ACGURYSWKMBDHVN"
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYBJOUXZ*"
OPERATIONS = ("sequence_lengths", "weighted_gc_fractions")
SUPPORTED_CASES = [
    ("sequence_lengths", "dna", DNA_ALPHABET),
    ("sequence_lengths", "rna", RNA_ALPHABET),
    ("sequence_lengths", "protein", PROTEIN_ALPHABET),
    ("weighted_gc_fractions", "dna", DNA_ALPHABET),
    ("weighted_gc_fractions", "rna", RNA_ALPHABET),
]

# Independently declared biological possibilities, not numeric GC weights.
_DNA_BASE_SETS = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

# Hand-derived fractions explicitly exercise every declared nucleotide symbol.
_SYMBOL_GC_FRACTIONS = [
    ("A", 0.0),
    ("C", 1.0),
    ("G", 1.0),
    ("T", 0.0),
    ("R", 1 / 2),
    ("Y", 1 / 2),
    ("S", 1.0),
    ("W", 0.0),
    ("K", 1 / 2),
    ("M", 1 / 2),
    ("B", 2 / 3),
    ("D", 1 / 3),
    ("H", 1 / 3),
    ("V", 2 / 3),
    ("N", 1 / 2),
]


def _rational_gc_oracle(value: str, kind: str) -> float:
    """Return exact set-derived GC probability rounded once to a float."""
    if not value:
        return 0.0
    base_sets = _DNA_BASE_SETS
    if kind == "rna":
        base_sets = {
            symbol.replace("T", "U"): frozenset(
                base.replace("T", "U") for base in bases
            )
            for symbol, bases in _DNA_BASE_SETS.items()
        }
    probabilities = (
        Fraction(len(base_sets[symbol] & {"G", "C"}), len(base_sets[symbol]))
        for symbol in value.upper()
    )
    return float(sum(probabilities, Fraction(0)) / len(value))


def _assert_gc_outputs(
    actual: list[float | None], expected: Sequence[float | None]
) -> None:
    """Compare all rows, preserving nulls and checking range and scalar types."""
    assert isinstance(actual, list)
    assert len(actual) == len(expected)
    for observed, wanted in zip(actual, expected, strict=True):
        if wanted is None:
            assert observed is None
        else:
            assert type(observed) is float
            assert math.isfinite(observed)
            assert 0.0 <= observed <= 1.0
            assert observed == pytest.approx(wanted, rel=1e-12, abs=1e-12)


def _deterministic_corpus(kind: str) -> list[str | None]:
    """Return exhaustive short words plus seeded, long and mixed-case rows."""
    alphabet = DNA_ALPHABET if kind == "dna" else RNA_ALPHABET
    values: list[str | None] = [
        "".join(symbols)
        for length in range(4)
        for symbols in product(alphabet, repeat=length)
    ]
    rng = random.Random(20261003)  # noqa: S311 - reproducible scientific fixture
    for length in [0, 1, 2, 3, 7, 15, 16, 31, 32, 63, 64, 255, 256, 4096]:
        for _ in range(6):
            values.append(
                "".join(rng.choice(alphabet + alphabet.lower()) for _ in range(length))
            )
    values.extend(symbol * 16384 for symbol in "bdhv")
    values[1:1] = [None, "", alphabet.lower(), None, alphabet.lower()]
    return values


@pytest.mark.parametrize("operation", OPERATIONS)
def test_public_metrics_are_compiled_native_callables(operation: str) -> None:
    """Reject Python replacements or a disguised native fallback module."""
    native = importlib.import_module("biov._native")
    origin = native.__file__
    assert origin is not None
    assert any(origin.endswith(suffix) for suffix in EXTENSION_SUFFIXES), origin
    assert Path(origin).parent == Path(biov.__file__).parent
    function = getattr(native, operation)
    assert callable(function)
    assert inspect.isbuiltin(function)
    assert getattr(biov, operation) is function
    assert operation in biov.__all__


@pytest.mark.parametrize(
    ("kind", "alphabet", "expected_length"),
    [
        ("dna", DNA_ALPHABET, 15),
        ("rna", RNA_ALPHABET, 15),
        ("protein", PROTEIN_ALPHABET, 27),
    ],
)
def test_lengths_count_all_validated_symbols(
    kind: str, alphabet: str, expected_length: int
) -> None:
    """Retain full lengths, row order, duplicates, empty strings and nulls."""
    values = [alphabet.lower(), None, "", alphabet, alphabet.lower(), None]
    before = values.copy()
    result = biov.sequence_lengths(values, kind=kind)
    assert result == [expected_length, None, 0, expected_length, expected_length, None]
    assert isinstance(result, list)
    assert result is not values
    assert all(value is None or type(value) is int for value in result)
    assert values == before
    result[0] = 0
    assert values == before


def test_protein_stops_each_count_as_one_symbol() -> None:
    """Count stored stop markers without terminating or translating proteins."""
    assert biov.sequence_lengths(
        ["*", "M*", "***", "B*Z*", "uoXj"], kind="protein"
    ) == [
        1,
        2,
        3,
        4,
        4,
    ]


@pytest.mark.parametrize(("symbol", "expected"), _SYMBOL_GC_FRACTIONS)
@pytest.mark.parametrize("kind", ["dna", "rna"])
def test_every_iupac_symbol_has_its_hand_derived_gc_weight(
    symbol: str, expected: float, kind: str
) -> None:
    """Check all canonical and ambiguous bases, including ASCII lowercase."""
    if kind == "rna":
        symbol = symbol.replace("T", "U")
    _assert_gc_outputs(
        biov.weighted_gc_fractions([symbol, symbol.lower()], kind=kind),
        [expected, expected],
    )


@pytest.mark.parametrize("kind", ["dna", "rna"])
def test_complete_gc_fixtures_preserve_all_positions_and_rows(kind: str) -> None:
    """Keep zero-GC and ambiguous symbols in each nonempty denominator."""
    values = ["GCN", "", None, "GDVV", "AWN", "SW", "NNNN", "acgtryswkmbdhvn"]
    if kind == "rna":
        values = [
            None if value is None else value.replace("t", "u") for value in values
        ]
    values.extend(["GCN", None])
    before = values.copy()
    result = biov.weighted_gc_fractions(values, kind=kind)
    _assert_gc_outputs(
        result, [5 / 6, 0.0, None, 2 / 3, 1 / 6, 1 / 2, 1 / 2, 1 / 2, 5 / 6, None]
    )
    assert result is not values
    assert values == before
    result[0] = 0.0
    assert values == before


@pytest.mark.parametrize("kind", ["dna", "rna"])
def test_full_gc_corpus_agrees_with_independent_rational_set_oracle(kind: str) -> None:
    """Compare complete outputs with probabilities derived from biological sets."""
    values = _deterministic_corpus(kind)
    before = values.copy()
    expected = [
        None if value is None else _rational_gc_oracle(value, kind) for value in values
    ]
    _assert_gc_outputs(biov.weighted_gc_fractions(values, kind=kind), expected)
    assert values == before


@pytest.mark.parametrize("kind", ["dna", "rna"])
def test_full_metric_corpus_agrees_with_pinned_biopython(
    kind: str, record_property: Callable[[str, object], None]
) -> None:
    """Use a separately implemented, explicitly versioned differential baseline."""
    record_property("biopython_reference_version", Bio.__version__)
    assert Bio.__version__ == REFERENCE_BIOPYTHON_VERSION, (
        "The scientific reference changed: install Biopython "
        f"{REFERENCE_BIOPYTHON_VERSION} from uv.lock, or review and update "
        "the reference version and sequence-contract provenance."
    )
    values = _deterministic_corpus(kind)
    expected_gc = [
        None if value is None else Bio.SeqUtils.gc_fraction(value, ambiguous="weighted")
        for value in values
    ]
    expected_lengths = [
        None if value is None else len(ReferenceSeq(value)) for value in values
    ]
    _assert_gc_outputs(biov.weighted_gc_fractions(values, kind=kind), expected_gc)
    assert biov.sequence_lengths(values, kind=kind) == expected_lengths


@pytest.mark.parametrize("kind", ["dna", "rna"])
def test_metrics_are_invariant_under_reverse_complement_and_partitioning(
    kind: str,
) -> None:
    """Check every result across transformations and native batch boundaries."""
    values = _deterministic_corpus(kind)
    before = values.copy()
    complements = biov.reverse_complements(values, kind=kind)
    lengths = biov.sequence_lengths(values, kind=kind)
    gc = biov.weighted_gc_fractions(values, kind=kind)
    assert biov.sequence_lengths(complements, kind=kind) == lengths
    _assert_gc_outputs(biov.weighted_gc_fractions(complements, kind=kind), gc)
    for split in [0, 1, 19, len(values)]:
        assert lengths == (
            biov.sequence_lengths(values[:split], kind=kind)
            + biov.sequence_lengths(values[split:], kind=kind)
        )
        _assert_gc_outputs(
            biov.weighted_gc_fractions(values[:split], kind=kind)
            + biov.weighted_gc_fractions(values[split:], kind=kind),
            gc,
        )
    assert values == before


@pytest.mark.parametrize(("operation", "kind", "alphabet"), SUPPORTED_CASES)
@pytest.mark.parametrize("values", [[], [None], [None, None], [""], [None, "", None]])
def test_metrics_keep_empty_and_null_distinct(
    operation: str, kind: str, alphabet: str, values: list[str | None]
) -> None:
    """Return a fresh list even when no nonempty sequence is present."""
    result = getattr(biov, operation)(values, kind=kind)
    expected = [None if value is None else 0 for value in values]
    assert result == expected
    assert result is not values
    if operation == "weighted_gc_fractions":
        _assert_gc_outputs(result, expected)
    else:
        assert all(value is None or type(value) is int for value in result)


@pytest.mark.parametrize(
    "values", [[], [None], [None, None], [""], ["ACD"], [None, "", None]]
)
def test_gc_rejects_protein_even_when_batch_has_no_sequence(
    values: list[str | None],
) -> None:
    """Reject unsupported biological kinds before relying on sequence content."""
    before = values.copy()
    with pytest.raises(TypeError):
        biov.weighted_gc_fractions(values, kind="protein")
    assert values == before


@pytest.mark.parametrize("operation", OPERATIONS)
@pytest.mark.parametrize("kind", ["", "DNA", "RNA", "Protein", " dna", "rna ", "other"])
@pytest.mark.parametrize("values", [[], [None], ["AC"]])
def test_unknown_kinds_do_not_get_inferred_or_repaired(
    operation: str, kind: str, values: list[str | None]
) -> None:
    """Reject unknown string kinds even for empty or all-null batches."""
    with pytest.raises(ValueError):
        getattr(biov, operation)(values, kind=kind)


@pytest.mark.parametrize("operation", OPERATIONS)
@pytest.mark.parametrize("kind", [None, 1, True, b"dna", ["dna"]])
def test_non_string_kinds_are_type_errors(operation: str, kind: object) -> None:
    """Distinguish invalid kind types from unknown string kinds."""
    with pytest.raises(TypeError):
        getattr(biov, operation)([], kind=kind)


@pytest.mark.parametrize("operation", OPERATIONS)
def test_metric_kind_is_required_and_keyword_only(operation: str) -> None:
    """Preserve the explicit public batch calling convention."""
    function = getattr(biov, operation)
    with pytest.raises(TypeError):
        function([])
    with pytest.raises(TypeError):
        function([], "dna")


@pytest.mark.parametrize("operation", OPERATIONS)
@pytest.mark.parametrize("values", [None, "ACGT", b"ACGT", ("ACGT",), {"ACGT": 1}])
def test_non_list_containers_are_rejected(operation: str, values: object) -> None:
    """Do not reinterpret a scalar sequence as a batch of symbols."""
    with pytest.raises(TypeError):
        getattr(biov, operation)(values, kind="dna")


@pytest.mark.parametrize("operation", OPERATIONS)
def test_iterators_and_dataframe_containers_are_not_silently_adapted(
    operation: str,
) -> None:
    """Keep list API input boundaries distinct from the pandas adapter."""
    for values in [
        iter(["AC", None]),
        (value for value in ["AC", None]),
        np.array(["AC", None]),
        pd.Series(["AC", None]),
        pd.array(["AC", None], dtype="biov.dna"),
    ]:
        with pytest.raises(TypeError):
            getattr(biov, operation)(values, kind="dna")


@pytest.mark.parametrize("operation", OPERATIONS)
@pytest.mark.parametrize(
    "invalid",
    [
        42,
        0.5,
        True,
        b"AC",
        ["AC"],
        {"sequence": "AC"},
        ReferenceSeq("AC"),
        float("nan"),
        pd.NA,
        pd.NaT,
    ],
)
def test_bad_element_types_reject_without_mutation(
    operation: str, invalid: object
) -> None:
    """Reject a bad row atomically, without stringifying or converting to null."""
    values = ["ac", None, invalid, "gt"]
    before = values.copy()
    with pytest.raises(TypeError):
        getattr(biov, operation)(values, kind="dna")
    assert all(
        actual is original for actual, original in zip(values, before, strict=True)
    )


@pytest.mark.parametrize(("operation", "kind", "alphabet"), SUPPORTED_CASES)
def test_exact_ascii_alphabet_is_validated_before_computing(
    operation: str, kind: str, alphabet: str
) -> None:
    """Test every ASCII symbol, including gaps, whitespace, digits and NUL."""
    accepted = set(alphabet + alphabet.lower())
    function = getattr(biov, operation)
    for code_point in range(128):
        symbol = chr(code_point)
        values = ["ac", None, f"A{symbol}C", ""]
        before = values.copy()
        if symbol in accepted:
            if operation == "sequence_lengths":
                assert function(values, kind=kind) == [2, None, 3, 0], repr(symbol)
            else:
                _assert_gc_outputs(
                    function(values, kind=kind),
                    [1 / 2, None, _rational_gc_oracle(f"A{symbol}C", kind), 0.0],
                )
        else:
            with pytest.raises(biov.SequenceValidationError):
                function(values, kind=kind)
        assert values == before


@pytest.mark.parametrize(("operation", "kind", "alphabet"), SUPPORTED_CASES)
@pytest.mark.parametrize(
    "code_points",
    [
        (0xE9,),
        (0xFF21,),
        (0xDF,),
        (0x17F,),
        (0x131,),
        (0xA0,),
        (0x1F9EC,),
        (0xD800,),
        (0xDFFF,),
        (0xD800, 0xDFFF),
    ],
)
def test_unicode_and_surrogates_raise_shared_sequence_error(
    operation: str, kind: str, alphabet: str, code_points: tuple[int, ...]
) -> None:
    """Never accept Unicode case folding or leak a surrogate encoding error."""
    invalid = "".join(chr(code_point) for code_point in code_points)
    values = ["ac", None, f"A{invalid}C", ""]
    before = values.copy()
    with pytest.raises(biov.SequenceValidationError):
        getattr(biov, operation)(values, kind=kind)
    assert values == before


@pytest.mark.parametrize(
    ("accessor", "operation", "kind", "alphabet", "dtype"),
    [
        ("length", "sequence_lengths", "dna", DNA_ALPHABET, "Int64"),
        ("length", "sequence_lengths", "rna", RNA_ALPHABET, "Int64"),
        ("length", "sequence_lengths", "protein", PROTEIN_ALPHABET, "Int64"),
        ("gc_fraction", "weighted_gc_fractions", "dna", DNA_ALPHABET, "Float64"),
        ("gc_fraction", "weighted_gc_fractions", "rna", RNA_ALPHABET, "Float64"),
    ],
)
@pytest.mark.parametrize("empty", [False, True])
def test_pandas_metrics_delegate_once_and_preserve_series_metadata(
    monkeypatch: pytest.MonkeyPatch,
    accessor: str,
    operation: str,
    kind: str,
    alphabet: str,
    dtype: str,
    empty: bool,
) -> None:
    """Adapt a complete pandas batch once, preserving index, name and nulls."""
    sequence_module = importlib.import_module("biov.seq")
    native = importlib.import_module("biov._native")
    native_function = getattr(native, operation)
    assert getattr(sequence_module, operation) is native_function
    inputs = [] if empty else [alphabet.lower(), None, "", alphabet.lower()]
    index = pd.Index([] if empty else [7, 3, 3, 9], name="row")
    source = pd.Series(inputs, index=index, name="sequences", dtype=f"biov.{kind}")
    before = source.copy(deep=True)
    adapted = [] if empty else [alphabet, None, "", alphabet]
    calls: list[tuple[list[str | None], str]] = []

    def record_native_call(values: list[str | None], *, kind: str) -> list[object]:
        calls.append((values.copy(), kind))
        return native_function(values, kind=kind)

    monkeypatch.setattr(sequence_module, operation, record_native_call)
    result = source.seq.length if accessor == "length" else source.seq.gc_fraction()
    assert calls == [(adapted, kind)]
    expected = pd.Series(
        native_function(adapted, kind=kind), index=index, dtype=dtype, name="sequences"
    )
    pd.testing.assert_series_equal(result, expected)
    pd.testing.assert_series_equal(source, before)


def test_native_and_pandas_gc_do_not_call_biopython(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep weighted GC computation out of the previous Python backend."""
    sequence_module = importlib.import_module("biov.seq")

    def forbidden_reference_call(*args: object, **kwargs: object) -> None:
        pytest.fail("Native or pandas GC computation called Biopython")

    monkeypatch.setattr(Bio.SeqUtils, "gc_fraction", forbidden_reference_call)
    if hasattr(sequence_module, "gc_fraction"):
        monkeypatch.setattr(sequence_module, "gc_fraction", forbidden_reference_call)
    for kind in ["dna", "rna"]:
        values = ["GCN", "", None, "BDHV", "GCN"]
        _assert_gc_outputs(
            biov.weighted_gc_fractions(values, kind=kind),
            [5 / 6, 0.0, None, 1 / 2, 5 / 6],
        )
        result = pd.Series(values, dtype=f"biov.{kind}").seq.gc_fraction()
        expected = pd.Series([5 / 6, 0.0, None, 1 / 2, 5 / 6], dtype="Float64")
        pd.testing.assert_series_equal(result, expected)


@pytest.mark.parametrize("values", [[], [None], [None, None], [""], ["ACD"]])
def test_pandas_gc_rejects_protein_without_data(values: list[str | None]) -> None:
    """Do not let pandas empty-batch adaptation evade the native kind contract."""
    source = pd.Series(values, dtype="biov.protein")
    before = source.copy(deep=True)
    with pytest.raises(TypeError):
        source.seq.gc_fraction()
    pd.testing.assert_series_equal(source, before)
