"""Check the BWA boundary with SAM fixtures and real BAM processing."""

import subprocess
from pathlib import Path

import pysam
import pytest

from biov import align_paired_reads, alignment


@pytest.mark.parametrize("name_sorted", [False, True])
def test_align_paired_reads_sorts_and_indexes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name_sorted: bool
) -> None:
    """Retain native options and use pysam to sort and optionally index SAM."""
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGTACGTAC\n")
    first, second = tmp_path / "r1.fq", tmp_path / "r2.fq"
    bam = tmp_path / "reads.bam"
    calls: list[tuple[str, ...]] = []

    def run_software(
        tool: str, arguments: tuple[str, ...], *, cwd: Path
    ) -> subprocess.CompletedProcess[bytes]:
        assert tool == "bwa"
        assert cwd == Path.cwd()
        calls.append(arguments)
        if arguments[0] == "mem":
            assert arguments[:3] == ("mem", "-t", "4")
            assert arguments[-3:] == (str(reference), str(first), str(second))
            Path(arguments[arguments.index("-o") + 1]).write_text(
                "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:10\n"
                "r2\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
                "r1\t0\tchr1\t2\t60\t4M\t*\t0\t0\tCGTA\tIIII\n"
            )
        return subprocess.CompletedProcess([tool, *arguments], 0)

    monkeypatch.setattr(alignment.settings, "execution_host", None)
    monkeypatch.setattr(alignment, "run_software", run_software)
    align_paired_reads(
        reference,
        first,
        second,
        bam,
        mem_args=("-t", "4"),
        sort_args=("-n", "-@", "4", "-m", "2G") if name_sorted else (),
        build_index=not name_sorted,
    )
    assert calls[0] == ("index", str(reference))
    assert len(calls) == 2
    with pysam.AlignmentFile(bam, "rb") as alignments:
        assert [read.query_name for read in alignments] == (
            ["r1", "r2"] if name_sorted else ["r2", "r1"]
        )
        assert alignments.has_index() is not name_sorted


@pytest.mark.parametrize("failed_step", ["index", "mem"])
def test_align_paired_reads_preserves_command_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failed_step: str
) -> None:
    """A failed BWA command cannot create a successful BAM artifact."""
    calls: list[str] = []

    def run_software(
        tool: str, arguments: tuple[str, ...], *, cwd: Path
    ) -> subprocess.CompletedProcess[bytes]:
        calls.append(arguments[0])
        return subprocess.CompletedProcess(
            [tool, *arguments], 7 if arguments[0] == failed_step else 0
        )

    monkeypatch.setattr(alignment.settings, "execution_host", None)
    monkeypatch.setattr(alignment, "run_software", run_software)
    bam = tmp_path / "reads.bam"
    with pytest.raises(subprocess.CalledProcessError):
        align_paired_reads(tmp_path / "ref.fa", tmp_path / "r1", tmp_path / "r2", bam)
    assert calls == (["index"] if failed_step == "index" else ["index", "mem"])
    assert not bam.exists()


def test_align_paired_reads_rejects_remote_temporary_paths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Do not pass local inputs and temporary outputs to SSH command dispatch."""
    monkeypatch.setattr(alignment.settings, "execution_host", "compute")
    with pytest.raises(ValueError, match="complete Python script"):
        align_paired_reads(Path("ref.fa"), Path("r1"), Path("r2"), Path("out.bam"))
