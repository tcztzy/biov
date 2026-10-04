"""GOATOOLS declaration is separate from installation and actual execution."""

import tomllib
from pathlib import Path

import pytest

from biov import environments


def test_goatools_is_an_independent_pinned_environment():
    """Declaring an environment neither installs it nor proves a scientific run."""
    with environments.manifest_source().open("rb") as source:
        pixi = tomllib.load(source)["tool"]["pixi"]
    feature = pixi["feature"]["goatools"]
    assert pixi["environments"]["goatools"] == {
        "features": ["goatools"],
        "no-default-feature": True,
    }
    assert feature["pypi-dependencies"] == {
        "goatools": "==1.6.5",
        "statsmodels": "==0.14.6",
    }
    assert feature["tasks"]["goatools"] == 'cd "$INIT_CWD" && goatools'
    assert "goatools" in environments.declared_environments()
    assert "goatools" not in pixi["feature"]["python"]["pypi-dependencies"]
    assert (
        Path(__file__).with_name("fixtures").joinpath("goatools", "tiny.obo").is_file()
    )


def test_saved_cli_output_matches_independent_exact_oracle():
    """Check saved real upstream output; this is not a new live-install claim."""
    import csv
    import importlib.util
    from fractions import Fraction

    script = Path(__file__).parents[1] / "scripts" / "validate_goatools.py"
    spec = importlib.util.spec_from_file_location("goatools_validation", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.fisher_two_sided(4, 0, 0, 6) == Fraction(1, 210)
    assert module.fisher_two_sided(0, 4, 6, 0) == Fraction(1, 210)
    assert module.fisher_two_sided(4, 0, 6, 0) == 1
    with (
        Path(__file__)
        .with_name("fixtures")
        .joinpath("goatools", "expected.tsv")
        .open() as stream
    ):
        reader = csv.DictReader(stream, delimiter="\t")
        rows = {row["# GO"]: row for row in reader}
    assert len(rows) == 3
    for go in ("GO:9000001", "GO:9000002"):
        assert float(rows[go]["p_uncorrected"]) == pytest.approx(1 / 210)
        assert float(rows[go]["p_bonferroni"]) == pytest.approx(1 / 70)
        assert float(rows[go]["p_fdr_bh"]) == pytest.approx(1 / 140)

    expected_ratios = {
        "GO:0008150": ("4/4", "10/10", "1, 2, 3, 4"),
        "GO:9000001": ("4/4", "4/10", "1, 2, 3, 4"),
        "GO:9000002": ("0/4", "6/10", ""),
    }
    for go, (study, population, members) in expected_ratios.items():
        assert rows[go]["ratio_in_study"] == study
        assert rows[go]["ratio_in_pop"] == population
        assert rows[go]["study_items"] == members
    root = rows["GO:0008150"]
    assert root["study_count"] == "4"
    for field in ("p_uncorrected", "p_bonferroni", "p_fdr_bh"):
        assert float(root[field]) == pytest.approx(1)


def test_validation_rejects_duplicate_go_rows(monkeypatch, tmp_path):
    """Reject repeated GO terms rather than silently overwriting a TSV row."""
    import importlib.util
    import subprocess

    script = Path(__file__).parents[1] / "scripts" / "validate_goatools.py"
    spec = importlib.util.spec_from_file_location("goatools_validation", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    fixture = Path(__file__).with_name("fixtures") / "goatools" / "expected.tsv"
    lines = fixture.read_text().splitlines(keepends=True)

    def duplicate_output(command, **kwargs):
        tmp_path.joinpath("results.tsv").write_text("".join([*lines, lines[1]]))
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(module.subprocess, "run", duplicate_output)
    monkeypatch.setattr(module.sys, "argv", [str(script), str(tmp_path)])
    with pytest.raises(ValueError, match="Duplicate GO row: GO:9000001"):
        module.main()
