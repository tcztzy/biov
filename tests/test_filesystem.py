"""Acceptance tests for identifier-backed fsspec reads."""

import json
import os
import subprocess
import sys
import tomllib
from pathlib import Path

import fsspec
import pytest
from test_artifacts import (
    ACCESSION,
    FASTA,
    UNIPROT_FASTA,
    UNIPROT_JSON,
    FakeDatasetsCli,
    FakeUniProtApi,
    _genome_package,
)

import biov.artifacts as artifacts_module
from biov import path, read_fasta, read_gff3

UNIPROT_URI = "uniprot://P42212"
REFSEQ_URI = f"refseq.gcf://{ACCESSION}"
GFF = b"##gff-version 3\nNC_003197.2\tRefSeq\tgene\t2\t5\t.\t+\t.\tID=gene1\n"


@pytest.mark.parametrize("configured", [False, True])
def test_cache_paths_respect_environment_and_fsspec_configuration(
    configured: bool, tmp_path: Path
) -> None:
    """Use the selected cache and preserve explicit native fsspec settings."""
    home = tmp_path / "home"
    expected = tmp_path / "files" if configured else home
    script = """
import os
import sys
from pathlib import Path
import fsspec
from biov.config import settings

assert settings.home == Path(os.environ["BIOV_HOME"])
memory = fsspec.filesystem("memory")
memory.pipe("/record.txt", b"cached data")
cache = fsspec.filesystem("filecache", fs=memory)
with cache.open("/record.txt", "rb") as stream:
    assert stream.read() == b"cached data"
    assert Path(stream.name).parent == Path(sys.argv[1])
if "FSSPEC_FILECACHE" in os.environ:
    assert Path(stream.name).name == "record.txt"
"""
    env = {
        k: v for k, v in os.environ.items() if not k.startswith(("FSSPEC_", "BIOV_"))
    }
    env["BIOV_HOME"] = str(home)
    env["FSSPEC_CONFIG_DIR"] = str(tmp_path / "fsspec")
    if configured:
        env["FSSPEC_FILECACHE"] = json.dumps(
            {"cache_storage": str(expected), "same_names": True}
        )
    subprocess.run(
        [sys.executable, "-c", script, str(expected)],
        env=env,
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )


def test_declared_protocols_match_manifest():
    """Register exactly the fsspec protocols the manifest declares."""
    project = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text())
    assert set(project["project"]["entry-points"]["fsspec.specs"]) == set(
        artifacts_module.artifact_capabilities()["namespaces"]
    )


@pytest.fixture
def providers(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> tuple[FakeUniProtApi, FakeDatasetsCli]:
    """Use the real artifact cache with deterministic provider downloads.

    Returns:
        Paired UniProt and NCBI Datasets fakes that record every download.
    """
    uniprot = FakeUniProtApi()
    datasets = FakeDatasetsCli(
        _genome_package(
            extra_members={f"ncbi_dataset/data/{ACCESSION}/genomic.gff": GFF}
        )
    )
    monkeypatch.setattr(artifacts_module.settings, "home", tmp_path)
    monkeypatch.setattr(
        artifacts_module, "download_uniprot_entry", uniprot.download_entry
    )
    monkeypatch.setattr(
        artifacts_module, "download_genome_package", datasets.download_genome_package
    )
    return uniprot, datasets


def test_fsspec_discovers_installed_protocols_without_importing_biov(
    providers: tuple[FakeUniProtApi, FakeDatasetsCli], tmp_path: Path
) -> None:
    """Discover all protocols and read cached files in a fresh process."""
    path(UNIPROT_URI)
    path(REFSEQ_URI)
    script = """
import json
import sys
from importlib.resources import files
import fsspec

assert "biov" not in sys.modules
from fsspec.registry import known_implementations

manifest = json.loads(files("biov").joinpath("assets/artifact_capabilities.json").read_text())
assert {name for name, spec in known_implementations.items()
        if spec["class"] == "biov.filesystem.BioVFileSystem"} == set(manifest["namespaces"])
for protocol in manifest["namespaces"]:
    assert known_implementations[protocol]["class"] == "biov.filesystem.BioVFileSystem"
contents = []
for uri in sys.argv[1:]:
    with fsspec.open(uri, "rt") as stream:
        contents.append(stream.read())
print(json.dumps(contents))
"""
    result = subprocess.run(
        [sys.executable, "-c", script, UNIPROT_URI, REFSEQ_URI],
        env={
            **{
                k: v
                for k, v in os.environ.items()
                if not k.startswith(("FSSPEC_", "BIOV_"))
            },
            "BIOV_HOME": str(tmp_path),
        },
        capture_output=True,
        text=True,
        check=True,
        timeout=30,
    )
    assert json.loads(result.stdout) == [UNIPROT_FASTA.decode(), FASTA.decode()]


@pytest.mark.parametrize(
    ("uri", "content", "sequence_id", "sequence"),
    [
        (
            UNIPROT_URI,
            UNIPROT_FASTA,
            "sp|P42212|GFP_AEQVI",
            "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDAT",
        ),
        (REFSEQ_URI, FASTA, "NC_003197.2", "ACGTGC"),
    ],
)
def test_fsspec_and_fasta_reader_reuse_the_artifact_cache(
    providers: tuple[FakeUniProtApi, FakeDatasetsCli],
    uri: str,
    content: bytes,
    sequence_id: str,
    sequence: str,
) -> None:
    """Keep text, binary, ranged and sequence reads on the same cached file."""
    with fsspec.open(uri, "rt") as stream:
        assert not isinstance(stream, list)
        assert stream.read() == content.decode()
    fs, fs_path = fsspec.core.url_to_fs(uri)
    assert fsspec.open_local(uri) == str(path(uri).path)
    assert fs.info(fs_path) == {"name": uri, "size": len(content), "type": "file"}
    assert fs.cat_file(fs_path, start=2, end=9) == content[2:9]
    assert fs.cat_file(fs_path, start=-5) == content[-5:]
    with fs.open(fs_path, "rb") as stream:
        assert stream.seek(3) == 3
        assert stream.read(5) == content[3:8]
        assert Path(stream.name) == path(uri).path
    records = read_fasta(uri)
    assert list(records) == [sequence_id]
    assert str(records[sequence_id].seq) == sequence
    uniprot, datasets = providers
    assert uniprot.download_calls == (
        [("protein_fasta", "P42212")] if uri == UNIPROT_URI else []
    )
    assert datasets.download_calls == ([ACCESSION] if uri == REFSEQ_URI else [])


def test_fsspec_artifact_option_selects_json_and_gff(
    providers: tuple[FakeUniProtApi, FakeDatasetsCli],
) -> None:
    """Select non-default artifacts through filesystem storage options."""
    with fsspec.open(UNIPROT_URI, "rt", artifact="entry_json") as stream:
        assert not isinstance(stream, list)
        assert stream.read() == UNIPROT_JSON.decode()
    with fsspec.open(REFSEQ_URI, "rb", artifact="annotation_gff3") as stream:
        assert not isinstance(stream, list)
        assert stream.read() == GFF
    annotation = read_gff3(REFSEQ_URI, storage_options={"artifact": "annotation_gff3"})
    assert len(annotation) == 1
    assert tuple(
        annotation.iloc[0][field] for field in ("seqid", "start", "end", "ID")
    ) == ("NC_003197.2", 1, 5, "gene1")
    uniprot, datasets = providers
    assert uniprot.download_calls == [("entry_json", "P42212")]
    assert datasets.download_calls == [ACCESSION]


@pytest.mark.parametrize("mode", ["wb", "ab", "xb", "r+b"])
def test_fsspec_rejects_mutating_modes_before_download(
    providers: tuple[FakeUniProtApi, FakeDatasetsCli], mode: str
) -> None:
    """Do not download or modify an artifact for a write request."""
    fs = fsspec.filesystem("uniprot")
    with pytest.raises(ValueError, match="read-only"):
        fs.open(UNIPROT_URI, mode)
    uniprot, datasets = providers
    assert uniprot.download_calls == []
    assert datasets.download_calls == []


def test_fsspec_preserves_provider_errors(
    providers: tuple[FakeUniProtApi, FakeDatasetsCli], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Surface the artifact provider exception without an error-text result."""
    failure = RuntimeError("provider unavailable")

    def fail(accession: str, destination: Path) -> None:
        raise failure

    monkeypatch.setattr(artifacts_module, "download_uniprot_entry", fail)
    with pytest.raises(RuntimeError) as raised, fsspec.open(UNIPROT_URI, "rb"):
        pass
    assert raised.value is failure


@pytest.mark.parametrize(
    "error",
    [RuntimeError("provider unavailable"), KeyboardInterrupt()],
    ids=["provider-error", "interrupt"],
)
def test_fsspec_exists_propagates_provider_failures(
    providers: tuple[FakeUniProtApi, FakeDatasetsCli],
    monkeypatch: pytest.MonkeyPatch,
    error: BaseException,
) -> None:
    """Never report a provider failure as a missing artifact."""
    fs = fsspec.filesystem("uniprot")

    def fail(accession: str, destination: Path) -> None:
        raise error

    monkeypatch.setattr(artifacts_module, "download_uniprot_entry", fail)
    with pytest.raises(type(error)) as raised:
        fs.exists(UNIPROT_URI)
    assert raised.value is error


def test_fsspec_exists_reports_only_genuinely_absent_files(
    providers: tuple[FakeUniProtApi, FakeDatasetsCli], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Resolve available artifacts and keep a provider "absent" answer as False."""
    fs = fsspec.filesystem("uniprot")
    assert fs.exists(UNIPROT_URI) is True

    def missing(accession: str, destination: Path) -> None:
        raise FileNotFoundError(f"UniProt has no entry {accession!r}")

    monkeypatch.setattr(artifacts_module, "download_uniprot_entry", missing)
    assert fs.exists("uniprot://P12345") is False
    assert fsspec.filesystem("encode").exists("encode://ENCSR163RYW") is False


@pytest.mark.parametrize("method", ["mkdir", "makedirs", "rmdir"])
def test_fsspec_directory_mutations_fail_before_provider_io(
    providers: tuple[FakeUniProtApi, FakeDatasetsCli], method: str
) -> None:
    """Refuse directory writes instead of silently succeeding."""
    fs = fsspec.filesystem("uniprot")

    with pytest.raises(NotImplementedError, match="read-only"):
        getattr(fs, method)("uniprot://created/nested")

    uniprot, datasets = providers
    assert uniprot.download_calls == []
    assert datasets.download_calls == []


def test_fsspec_preserves_unsupported_artifact_error(
    providers: tuple[FakeUniProtApi, FakeDatasetsCli],
) -> None:
    """Keep unsupported namespace/kind errors distinct from provider failures."""
    with (
        pytest.raises(ValueError),
        fsspec.open(UNIPROT_URI, "rb", artifact="annotation_gff3"),
    ):
        pass
    uniprot, datasets = providers
    assert uniprot.download_calls == []
    assert datasets.download_calls == []


def test_cache_home_expands_tilde(monkeypatch):
    """Expand a leading tilde in the configured cache home."""
    from biov.config import Settings

    monkeypatch.setenv("BIOV_HOME", "~/biov-cache")
    assert Settings().home == Path("~/biov-cache").expanduser().resolve()
