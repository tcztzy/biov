"""Regression coverage for migrated sequence and annotation boundaries."""

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import crisprprimer
from biov import BioDataFrame
from crisprprimer import id_converter
from crisprprimer.nuclease import AsCas12, Nuclease, Nucleases, SpCas9


@pytest.mark.parametrize(
    "sequence,start,end,strand,nuclease,expected",
    [
        ("TT" + "ACGT" * 5 + "AGG", 2, 22, "+", SpCas9, "AGG"),
        ("CCT" + "ACGT" * 5, 3, 23, "-", SpCas9, "AGG"),
        ("TTTA" + "ACGT" * 5, 4, 24, "+", AsCas12, "TTTA"),
        ("ACGT" * 5 + "TAAA", 0, 20, "-", AsCas12, "TTTA"),
    ],
)
def test_native_alignment_pam_orientation(
    sequence: str, start: int, end: int, strand: str, nuclease: Nuclease, expected: str
) -> None:
    """Retain both PAM orientations for existing nuclease presets."""
    row = pd.Series({"tName": "ref", "tStart": start, "tEnd": end, "strand": strand})
    assert (
        crisprprimer.get_pam(row, {"ref": SeqRecord(Seq(sequence))}, nuclease)
        == expected
    )


@pytest.mark.parametrize("strand,expected", [("+", "AGTC" * 5), ("-", "GACT" * 5)])
def test_native_protospacer_orientation(strand: str, expected: str) -> None:
    """Preserve the native PSL context when slicing a sequence record."""
    row = pd.Series(
        {
            "tName": "ref",
            "tStart": 3,
            "tEnd": 23,
            "qStart": 0,
            "qEnd": 20,
            "blockCount": 1,
            "strand": strand,
        }
    )
    records = {"ref": SeqRecord(Seq("CCC" + "AGTC" * 5 + "GGG"))}
    assert crisprprimer.get_protospacer(row, records) == expected


def test_alignment_context_rejects_missing_sequence() -> None:
    """Report missing reference bases instead of consuming a null sequence."""
    row = pd.Series({"tName": "ref", "tStart": 0, "tEnd": 20, "strand": "+"})
    records = {"ref": SeqRecord(None)}
    with pytest.raises(ValueError, match="ref has no sequence"):
        crisprprimer.get_pam(row, records)
    with pytest.raises(ValueError, match="ref has no sequence"):
        crisprprimer.get_protospacer(row, records)


def test_region_lookup_handles_missing_cds_and_unsupported_rap() -> None:
    """Empty MSU results terminate safely and unimplemented RAP lookup is explicit."""
    annotation = BioDataFrame(columns=["type", "ID", "seqid", "start", "end"])
    assert (
        crisprprimer._crispr_for_one_region(
            "LOC_Os01g00001",
            annotation,
            Nucleases.SpCas9,
            (37, 70),
            {},
            "unused.fa",
            (),
            (("", ""), ("", "")),
        )
        is None
    )
    with pytest.raises(ValueError, match="RAP annotation lookup is not implemented"):
        crisprprimer._crispr_for_one_region(
            "Os01g0000100",
            annotation,
            Nucleases.SpCas9,
            (37, 70),
            {},
            "unused.fa",
            (),
            (("", ""), ("", "")),
        )


def test_public_identifier_conversion_retains_mapping(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Public RAP/MSU lookup does not depend on the unpublished NAU table."""
    monkeypatch.setattr(id_converter, "RAP2MSU", {"Os01g0000100": ["LOC_Os01g00001.1"]})
    monkeypatch.setattr(id_converter, "MSU2RAP", {"LOC_Os01g00001.1": "Os01g0000100"})
    monkeypatch.setattr(id_converter, "_mapping_loaded", True)
    assert id_converter.convert("Os01g0000100") == ["LOC_Os01g00001.1"]
    assert id_converter.convert("LOC_Os01g00001") == "Os01g0000100"
    assert not hasattr(id_converter, "NAU2MSU")
    assert not hasattr(id_converter, "MSU2NAU")


def test_unsupported_identifier_rejected_before_download(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unsupported identifier schemes cannot silently return an empty mapping."""

    def unexpected_download() -> None:
        pytest.fail("Unsupported identifier triggered a mapping download")

    monkeypatch.setattr(id_converter, "_load_mapping", unexpected_download)
    with pytest.raises(ValueError, match="not RAP id or MSU id"):
        id_converter.convert("unsupported-id")
