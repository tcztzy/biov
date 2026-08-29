"""Tests for scalar and typed pandas sequence APIs."""

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal
from pydantic import BaseModel, TypeAdapter

from biov.seq import Seq, SequenceValidationError


class SequenceModel(BaseModel):
    """Pydantic fixture containing BioV's scalar sequence type."""

    seq: Seq


def test_seq_validation_and_serialization() -> None:
    """Keep the existing Biopython/Pydantic scalar contract."""
    acgt = TypeAdapter(Seq).validate_python("ACGT")
    assert acgt == "ACGT"
    assert acgt.reverse_complement() == "ACGT"
    model = SequenceModel.model_validate({"seq": "ACGT"})
    assert model.seq == acgt
    assert model.model_dump() == {"seq": "ACGT"}
    assert model.model_dump_json() == '{"seq":"ACGT"}'


def test_sequence_dtype_is_explicit_nullable_and_normalized() -> None:
    """Require semantic dtype instead of guessing from strings."""
    dna = pd.Series(["acgt", None, ""], dtype="biov.dna", name="sequence")

    assert str(dna.dtype) == "biov.dna"
    assert dna.iloc[0] == "ACGT"
    assert dna.iloc[1] is pd.NA
    assert dna.iloc[2] == ""
    assert_series_equal(
        dna.seq.length,
        pd.Series([4, pd.NA, 0], dtype="Int64", name="sequence"),
    )
    with pytest.raises(AttributeError, match="BioV sequence dtype"):
        _ = pd.Series(["ACGT"]).seq


def test_sequence_array_supports_indexing_assignment_and_concat() -> None:
    """Provide the minimum pandas ExtensionArray behavior callers rely on."""
    series = pd.Series(["AC", None, "GT"], dtype="biov.dna")
    selected = series.iloc[[2, 0]]
    copied = series.copy()
    copied.iloc[0] = "nn"
    combined = pd.concat([selected, copied], ignore_index=True)

    assert selected.tolist() == ["GT", "AC"]
    assert copied.tolist() == ["NN", pd.NA, "GT"]
    assert str(combined.dtype) == "biov.dna"
    with pytest.raises(SequenceValidationError, match="invalid DNA symbol"):
        copied.iloc[0] = "AU"


def test_sequence_array_supports_nullable_comparison_and_string_cast() -> None:
    """Interoperate with standard pandas comparison and string conversion."""
    series = pd.Series(["AC", None, "GT"], dtype="biov.dna")

    assert_series_equal(
        series == "AC",
        pd.Series([True, pd.NA, False], dtype="boolean"),
    )
    assert_series_equal(
        series != "AC",
        pd.Series([False, pd.NA, True], dtype="boolean"),
    )
    assert_series_equal(
        series.astype("string"),
        pd.Series(["AC", pd.NA, "GT"], dtype="string"),
    )


@pytest.mark.parametrize(
    ("dtype", "value", "message"),
    [
        ("biov.dna", "ACGU", "invalid DNA symbol"),
        ("biov.rna", "ACGT", "invalid RNA symbol"),
        ("biov.protein", "ACD?", "invalid protein symbol"),
        ("biov.dna", 42, "must be a string or missing"),
    ],
)
def test_sequence_validation_errors_are_stable(
    dtype: str, value: object, message: str
) -> None:
    """Reject wrong scalar types and alphabet violations predictably."""
    with pytest.raises(SequenceValidationError, match=message):
        pd.Series([value], dtype=dtype)


def test_declared_iupac_symbols_are_storable() -> None:
    """Store extended IUPAC symbols without treating them as canonical."""
    assert pd.Series(["RYSWKMBDHVN"], dtype="biov.dna").iloc[0] == "RYSWKMBDHVN"
    assert pd.Series(["RYSWKMBDHVN"], dtype="biov.rna").iloc[0] == "RYSWKMBDHVN"
    assert pd.Series(["BXZJUO*"], dtype="biov.protein").iloc[0] == "BXZJUO*"


@pytest.mark.parametrize(
    ("dtype", "value", "expected"),
    [
        ("biov.dna", "ACGTRYN", "NRYACGT"),
        ("biov.rna", "ACGURYN", "NRYACGU"),
    ],
)
def test_reverse_complement_preserves_sequence_dtype_and_nulls(
    dtype: str, value: str, expected: str
) -> None:
    """Use Biopython complement rules for DNA and RNA IUPAC symbols."""
    result = pd.Series([value, None], dtype=dtype).seq.reverse_complement()

    assert result.tolist() == [expected, pd.NA]
    assert str(result.dtype) == dtype


def test_gc_fraction_uses_weighted_iupac_and_handles_empty_and_null() -> None:
    """Define GC fractions for ambiguous, empty, and missing sequences."""
    result = pd.Series(["GCN", "", None], dtype="biov.dna").seq.gc_fraction()

    assert result.iloc[0] == pytest.approx(5 / 6)
    assert result.iloc[1] == pytest.approx(0.0)
    assert result.iloc[2] is pd.NA
    assert str(result.dtype) == "Float64"


def test_nucleic_acid_translation_returns_typed_protein() -> None:
    """Translate DNA/RNA complete codons and preserve missing values."""
    dna = pd.Series(["ATGGCC", "ATGTAA", "ATGNNN", None], dtype="biov.dna")
    rna = pd.Series(["AUGGCC"], dtype="biov.rna")

    translated = dna.seq.translate()
    stopped = dna.seq.translate(to_stop=True)

    assert translated.tolist() == ["MA", "M*", "MX", pd.NA]
    assert stopped.tolist() == ["MA", "M", "MX", pd.NA]
    assert rna.seq.translate().tolist() == ["MA"]
    assert str(translated.dtype) == "biov.protein"
    with pytest.raises(SequenceValidationError, match="complete codons"):
        pd.Series(["ATGGC"], dtype="biov.dna").seq.translate()


def test_protein_analysis_and_composition_use_biopython() -> None:
    """Return nullable scalar analyses and stable 20-column percentages."""
    protein = pd.Series(["ACD", None], dtype="biov.protein", name="protein")

    mass = protein.seq.molecular_weight()
    pi = protein.seq.isoelectric_point()
    composition = protein.seq.amino_acid_composition()

    assert mass.iloc[0] == pytest.approx(307.3235)
    assert mass.iloc[1] is pd.NA
    assert pi.iloc[0] == pytest.approx(4.299494743347168)
    assert pi.iloc[1] is pd.NA
    expected = pd.DataFrame(
        [
            {
                amino_acid: (100 / 3 if amino_acid in "ACD" else 0.0)
                for amino_acid in "ACDEFGHIKLMNPQRSTVWY"
            },
            dict.fromkeys("ACDEFGHIKLMNPQRSTVWY", pd.NA),
        ],
        dtype="Float64",
    )
    assert_frame_equal(composition, expected)
    assert_frame_equal(
        pd.Series([], dtype="biov.protein").seq.amino_acid_composition(),
        pd.DataFrame(columns=pd.Index(list("ACDEFGHIKLMNPQRSTVWY")), dtype="Float64"),
    )


@pytest.mark.parametrize("value", ["X", "M*", ""])
@pytest.mark.parametrize(
    "method", ["molecular_weight", "isoelectric_point", "amino_acid_composition"]
)
def test_protein_analysis_rejects_extended_or_empty_sequences(
    value: str, method: str
) -> None:
    """Wrap backend limitations in one BioV-owned error contract."""
    accessor = pd.Series([value], dtype="biov.protein").seq

    with pytest.raises(SequenceValidationError, match="non-empty canonical"):
        getattr(accessor, method)()


def test_sequence_accessor_rejects_operations_for_the_wrong_kind() -> None:
    """Keep nucleotide and protein method families separate."""
    with pytest.raises(TypeError, match="DNA or RNA"):
        pd.Series(["ACD"], dtype="biov.protein").seq.reverse_complement()
    with pytest.raises(TypeError, match="protein"):
        pd.Series(["ACG"], dtype="biov.dna").seq.molecular_weight()
