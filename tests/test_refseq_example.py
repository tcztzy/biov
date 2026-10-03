"""Guard checks for the explicitly scoped, independent RefSeq example reader.

The real provider download is an opt-in documented experiment, not a network
fixture in CI. These tests exercise local failure modes without BioV imports.
"""

import runpy
import subprocess
import sys
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "inspect_refseq_example.py"
_READER = runpy.run_path(str(_SCRIPT))


@pytest.mark.parametrize(
    "name",
    ["", ".", "..", "../outside", "/absolute", "a/../b", "a//b", r"a\b", "C:drive"],
)
def test_refseq_example_rejects_unsafe_paths(tmp_path: Path, name: str) -> None:
    """Reject traversal and ambiguous platform path syntax before reading data."""
    with pytest.raises(ValueError, match="Unsafe relative path"):
        _READER["safe_file"](tmp_path, name)


def test_refseq_example_confines_links_and_preserves_fasta(tmp_path: Path) -> None:
    """Keep ordinary relative files usable and reject a link outside the root."""
    root = tmp_path / "package"
    root.mkdir()
    fasta = root / "example.fna"
    fasta.write_text(
        ">NC_000001.2 description\nACGT\nNN\n>second\nTT\n", encoding="utf-8"
    )
    assert _READER["safe_file"](root, "example.fna") == fasta
    assert _READER["read_fasta"](fasta) == {"NC_000001.2": "ACGTNN", "second": "TT"}
    outside = tmp_path / "outside.fna"
    outside.write_text(">outside\nA\n", encoding="utf-8")
    link = root / "escape.fna"
    link.symlink_to(outside)
    with pytest.raises(ValueError, match="Unsafe or missing file"):
        _READER["safe_file"](root, "escape.fna")


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        (">duplicate\nA\n>duplicate\nT\n", "Duplicate FASTA ID"),
        (">\nA\n", "Empty FASTA header"),
        ("A\n>later\nT\n", "Sequence precedes first"),
        ("", "Empty FASTA file"),
    ],
)
def test_refseq_example_rejects_ambiguous_fasta(
    tmp_path: Path, contents: str, message: str
) -> None:
    """Fail explicitly instead of silently losing duplicate or malformed records."""
    path = tmp_path / "invalid.fna"
    path.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        _READER["read_fasta"](path)


def test_refseq_example_checks_survive_python_optimization(tmp_path: Path) -> None:
    """The standalone reader's validation must remain active with python -O."""
    result = subprocess.run(
        [
            sys.executable,
            "-I",
            "-O",
            "-c",
            "import runpy, sys; runpy.run_path(sys.argv[1])['require'](False, 'guard-active')",
            str(_SCRIPT),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    assert result.returncode != 0
    assert "ValueError: guard-active" in result.stderr
