"""Regression checks for the Azimuth input and backend contract."""

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from biov import azimuth as azimuth_mod


def test_sequence_record_normalizes_sequence() -> None:
    """Sequence record normalizes sequence."""
    record = azimuth_mod.SequenceRecord(sequence=" acagctgatctccagatatgaccatgggtt ")
    assert record.sequence == "ACAGCTGATCTCCAGATATGACCATGGGTT"


def test_sequence_record_requires_complete_position_pair() -> None:
    """Sequence record requires complete position pair."""
    with pytest.raises(ValidationError):
        azimuth_mod.SequenceRecord(
            sequence="ACAGCTGATCTCCAGATATGACCATGGGTT",
            aa_cut=12,
        )


def test_load_records_from_csv_accepts_30mer_alias(tmp_path: Path) -> None:
    """Load records from csv accepts 30mer alias."""
    input_path = tmp_path / "guides.csv"
    input_path.write_text(
        "id,30mer,aa_cut,percent_peptide\n"
        "guide-1,ACAGCTGATCTCCAGATATGACCATGGGTT,12,0.25\n",
        encoding="utf-8",
    )
    records = azimuth_mod.load_records_from_path(input_path)
    assert len(records) == 1
    assert records[0].record_id == "guide-1"
    assert records[0].sequence == "ACAGCTGATCTCCAGATATGACCATGGGTT"
    assert records[0].aa_cut == 12
    assert records[0].percent_peptide == pytest.approx(0.25)


def test_load_records_from_jsonl(tmp_path: Path) -> None:
    """Load records from jsonl."""
    input_path = tmp_path / "guides.jsonl"
    input_path.write_text(
        json.dumps({"id": "guide-1", "sequence": "ACAGCTGATCTCCAGATATGACCATGGGTT"})
        + "\n",
        encoding="utf-8",
    )
    records = azimuth_mod.load_records_from_path(input_path)
    assert [record.record_id for record in records] == ["guide-1"]


def test_choose_backend_auto_prefers_docker(monkeypatch: pytest.MonkeyPatch) -> None:
    """Choose backend auto prefers docker."""
    monkeypatch.setattr(azimuth_mod, "docker_image_available", lambda image: True)
    assert azimuth_mod.choose_backend("auto", "azimuth:latest") == "docker"


def test_choose_backend_auto_falls_back_to_fork(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Choose backend auto falls back to fork."""
    monkeypatch.setattr(azimuth_mod, "docker_image_available", lambda image: False)
    assert azimuth_mod.choose_backend("auto", "azimuth:latest") == "fork"
