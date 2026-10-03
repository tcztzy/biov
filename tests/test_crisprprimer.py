"""Check the migrated package and its native BLAT execution boundary."""

import gzip
import importlib
import subprocess
import tomllib
from pathlib import Path

import fsspec.config
import pandas as pd
import pytest
from Bio.Seq import MutableSeq, Seq
from typer.testing import CliRunner

import crisprprimer
import crisprprimer.__main__ as cli
from crisprprimer.__main__ import app
from crisprprimer.docker import parse_html_report
from crisprprimer.id_converter import guess_id_system
from crisprprimer.nuclease import AsCas12, Nuclease, Nucleases, SpCas9
from crisprprimer.score import cfd_score, pair_score


def test_package_resources_and_entry_points(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Retain the existing imports, bundled scores, and CLI names in BioV."""
    assert guess_id_system("Os01g0183666") == "RAP"
    assert guess_id_system("LOC_Os12g16350.10") == "MSU"
    nuclease = Nuclease.model_validate(
        {"cut_sites": -3, "pam": "NGG", "spacer_range": (-20, 0)}
    )
    assert nuclease.prototype == "N" * 20 + "NGG"
    assert cfd_score("A" * 20, "A" * 20, "GG") == pytest.approx(1.0)
    assert crisprprimer.RESTRICTION_ENZYMES["BsaI"] == "GGTCTC"
    with pytest.raises(ValueError, match="CFD requires"):
        cfd_score("A", "A", "GG")
    assert CliRunner().invoke(app, ["--help"]).exit_code == 0
    project = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text())
    assert project["project"]["scripts"]["crisprprimer"] == "crisprprimer.__main__:run"
    assert (
        project["project"]["scripts"]["crisprprimer-docker"]
        == "crisprprimer.docker:main"
    )
    assert project["project"]["scripts"]["biov-azimuth"] == "biov.azimuth:main"
    selected_cache = str(tmp_path / "configured-cache")
    monkeypatch.setitem(
        fsspec.config.conf["filecache"], "cache_storage", selected_cache
    )
    importlib.reload(crisprprimer)
    assert fsspec.config.conf["filecache"]["cache_storage"] == selected_cache


@pytest.mark.parametrize("returncode,empty", [(0, False), (0, True), (7, False)])
def test_blat_uses_native_execution_and_preserves_failures(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, returncode: int, empty: bool
) -> None:
    """Preserve PSL fields, compressed input, cache reuse, and process failures."""
    reference = tmp_path / "reference.fa.gz"
    reference.write_bytes(gzip.compress(b">reference\nACGT\n"))
    query = "ACGT"
    calls: list[tuple[str, ...]] = []
    monkeypatch.setattr(crisprprimer.settings, "execution_host", None)

    def run_software(
        tool: str, arguments: tuple[str, ...], *, cwd: Path
    ) -> subprocess.CompletedProcess[bytes]:
        assert tool == "blat"
        assert cwd == Path.cwd()
        assert arguments[:5] == (
            "-noHead",
            "-minMatch=0",
            "-minScore=18",
            "-stepSize=5",
            "-fine",
        )
        assert Path(arguments[-3]).read_bytes() == b">reference\nACGT\n"
        assert Path(arguments[-2]).read_text() == ">ACGT\nACGT\n"
        calls.append(arguments)
        Path(arguments[-1]).write_text(
            ""
            if empty
            else "4\t0\t0\t0\t0\t0\t0\t0\t+\tACGT\t4\t0\t4\tref\t4\t0\t4\t1\t4,\t0,\t0,\n"
        )
        return subprocess.CompletedProcess([tool, *arguments], returncode)

    monkeypatch.setattr(crisprprimer, "run_software", run_software)
    cache = tmp_path / "cache"
    if returncode:
        with pytest.raises(subprocess.CalledProcessError):
            crisprprimer._blat(str(reference), query, cache)
        assert not list(cache.glob("**/*.parquet"))
        return
    result = crisprprimer._blat(str(reference), query, cache)
    assert list(result.columns) == crisprprimer._BLAT_COLUMNS
    assert result.empty is empty
    if not empty:
        assert result.iloc[0]["qName"] == query
        assert result.iloc[0]["blockSizes"] == "4,"
        assert (
            crisprprimer._blat(str(reference), query, cache).iloc[0]["qName"] == query
        )
        assert len(calls) == 1


def test_blat_rejects_remote_temporary_paths(monkeypatch: pytest.MonkeyPatch) -> None:
    """Do not send caller-local temporary files to a separately configured host."""
    monkeypatch.setattr(crisprprimer.settings, "execution_host", "compute")
    with pytest.raises(ValueError, match="complete Python script"):
        crisprprimer._blat("unused.fa", "ACGT")


def test_parse_legacy_html_report(tmp_path: Path) -> None:
    """Retain the existing Docker report's public column mapping."""
    report = tmp_path / "A.html"
    report.write_text(
        "<table><tr><th>rank</th><th>ID</th><th>CHROM</th>"
        "<th>PAM start</th><th>PAM end</th><th>PAM Out</th>"
        "<th>Spacer 1</th><th>Spacer 2</th><th>Score</th></tr>"
        "<tr><td>1</td><td>gene-1</td><td>Chr1</td><td>100</td><td>102</td>"
        "<td>+</td><td>ACGT</td><td>TGCA</td><td>87</td></tr></table>"
    )
    hits = parse_html_report(report)
    assert len(hits) == 1
    assert hits[0].chrom == "Chr1"
    assert hits[0].score == 87


@pytest.mark.parametrize("name,model", [("SpCas9", SpCas9), ("AsCas12a", AsCas12)])
def test_nuclease_names_preserve_models(name: str, model: Nuclease) -> None:
    """Retain the public names and original scientific configurations."""
    selected = Nucleases[name]
    assert selected.value == name
    assert Nucleases(name) is selected
    assert selected.model is model


@pytest.mark.parametrize(
    "system,expected",
    [
        (None, Nucleases.SpCas9),
        ("SpCas9", Nucleases.SpCas9),
        ("AsCas12a", Nucleases.AsCas12a),
    ],
)
def test_cli_passes_selected_nuclease(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    system: str | None,
    expected: Nucleases,
) -> None:
    """Accept default and explicit CLI names without starting an analysis."""
    calls: list[Nucleases] = []

    def design(
        *, regions: list[str] | None, system: Nucleases, blocks_dir: Path
    ) -> None:
        assert regions == ["Chr1:1..100"]
        assert blocks_dir == tmp_path / "output"
        calls.append(system)

    monkeypatch.setattr(cli, "crisprprimer", design)
    args = [str(tmp_path / "output"), "-r", "Chr1:1..100"]
    if system is not None:
        args.extend(["--system", system])
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output
    assert calls == [expected]


def test_cfd_retains_mismatch_and_pam_factors() -> None:
    """Keep the bundled positional mismatch product and unsupported-PAM zero."""
    spacer = "A" * 11 + "TT" + "A" * 7
    assert cfd_score(spacer, "A" * 20, "AG") == pytest.approx(
        0.8 * 0.692307692 * 0.259259259
    )
    assert cfd_score("A" * 20, "A" * 20, "NN") == pytest.approx(0.0)


@pytest.mark.parametrize(
    "strand,start,end,expected",
    [
        ("+", 110, 176, 0.785),
        ("-", 110, 176, 0.5),
        ("+", 100, 300, 0.5),
        ("+", 110, 112, 0.0),
    ],
)
def test_pair_score_preserves_distance_and_length(
    strand: str, start: int, end: int, expected: float
) -> None:
    """Retain strand-aware CDS distance, full coverage, and short-pair scores."""
    cds = pd.DataFrame({"start": [100], "end": [300], "strand": [strand]})
    row = pd.Series({"start": start, "end": end})
    assert pair_score(row, cds, dist=3) == pytest.approx(expected)


@pytest.mark.parametrize("sequence_type", [Seq, MutableSeq])
def test_nuclease_accepts_biopython_sequences(
    sequence_type: type[Seq] | type[MutableSeq],
) -> None:
    """Keep native immutable and mutable sequence inputs and forward coordinates."""
    sequence = sequence_type("A" * 20 + "AGG")
    found = SpCas9.find_spacers_on(sequence)
    assert list(found["start"]) == [0]
    assert list(found["end"]) == [20]
    assert list(found["pam"]) == ["AGG"]
    assert list(found["strand"]) == ["+"]
    assert list(found["spacer"]) == ["A" * 20]
    assert str(sequence) == "A" * 20 + "AGG"
