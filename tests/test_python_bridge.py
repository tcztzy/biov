"""Distribution/interpreter pairing checks without importing BioV's runtime."""

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest


@pytest.fixture
def bridge():
    """Load only the private bridge, independent of scientific dependencies.

    Returns:
        The loaded private bridge module.
    """
    path = Path(__file__).parents[1] / "src" / "biov" / "_bridge.py"
    spec = importlib.util.spec_from_file_location("bridge_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pairing_accepts_only_recorded_real_binary(bridge, tmp_path, monkeypatch):
    """A same-named package cannot authorize another binary in its sibling venv."""
    binary = tmp_path / "bin" / "biov"
    binary.parent.mkdir()
    binary.write_bytes(b"native binary fixture")
    entry = Path("../../../bin/biov")
    monkeypatch.setattr(
        bridge,
        "distribution",
        lambda _: SimpleNamespace(files=[entry], locate_file=lambda _: binary),
    )
    bridge.validate_pairing(binary)
    other = tmp_path / "other-biov"
    other.write_bytes(binary.read_bytes())
    with pytest.raises(ValueError, match="does not own"):
        bridge.validate_pairing(other)


def test_pairing_rejects_missing_record_and_absent_binary(
    bridge, tmp_path, monkeypatch
):
    """Metadata absence is an error, never permission to fall back to PATH."""
    binary = tmp_path / "biov"
    binary.write_bytes(b"native binary fixture")
    monkeypatch.setattr(bridge, "distribution", lambda _: SimpleNamespace(files=None))
    with pytest.raises(ValueError, match="does not own"):
        bridge.validate_pairing(binary)
    with pytest.raises(FileNotFoundError):
        bridge.validate_pairing(tmp_path / "missing")


def test_bridge_requires_caller_before_importing_legacy_cli(
    bridge, monkeypatch, capsys
):
    """Direct invocation without a recorded native caller exits cleanly."""
    monkeypatch.setattr(sys, "argv", ["biov._bridge"])
    monkeypatch.setattr(bridge, "validate_pairing", Mock(side_effect=AssertionError))
    with pytest.raises(SystemExit) as error:
        bridge.main()
    assert error.value.code == 2
    assert capsys.readouterr().out == ""
