"""Acceptance tests for identifier-backed environment-local artifacts."""

import io
import json
import os
import subprocess
import zipfile
from pathlib import Path
from typing import Any

import pytest

import biov.artifacts as artifacts_module
from biov.artifacts import (
    DATASETS_CLI_TIMEOUT_SECONDS,
    DATASETS_INCLUDE,
    ArtifactNotFoundError,
    ArtifactPackageError,
    ArtifactServiceError,
    UnsupportedArtifactError,
    download_genome_package,
    download_uniprot_entry,
    genome_summary,
    open,
    path,
)
from biov.identifiers import IdentifierSyntaxError, parse_identifier

ACCESSION = "GCF_000006945.2"
VERSIONLESS_ACCESSION = "GCF_000006945"
FASTA_NAME = f"{ACCESSION}_ASM694v2_genomic.fna"
FASTA = b">NC_003197.2 chromosome\nACGTGC\n"
UNIPROT_ACCESSION = "P42212"
UNIPROT_FASTA = (
    b">sp|P42212|GFP_AEQVI Green fluorescent protein OS=Aequorea victoria "
    b"OX=6100 GN=GFP PE=1 SV=1\nMSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDAT\n"
)
UNIPROT_JSON = (
    b'{\n  "entryType": "UniProtKB reviewed (Swiss-Prot)",\n'
    b'  "primaryAccession": "P42212",\n  "uniProtkbId": "GFP_AEQVI",\n'
    b'  "uniProtKBCrossReferences": [\n'
    b'    {"database": "PDB", "id": "1EMA", "properties": []},\n'
    b'    {"database": "PDB", "id": "1GFL", "properties": []}\n'
    b"  ]\n}\n"
)
GENOME_SUMMARY = (
    '{"reports":[{"accession":"GCF_000006945.2",'
    '"organism":{"organism_name":"Salmonella enterica"}}],"total_count":1}\n'
)


class FakeDatasetsCli:
    """Deterministic replacement for the external NCBI Datasets CLI."""

    def __init__(self, package: bytes | None = None) -> None:
        """Initialize a fake with one prepared data package."""
        self.package = package or _genome_package()
        self.download_calls: list[str] = []

    def download_genome_package(self, accession: str, destination: Path) -> None:
        """Write the prepared ZIP and record the requested accession."""
        self.download_calls.append(accession)
        destination.write_bytes(self.package)


class FakeUniProtApi:
    """Deterministic replacement for the UniProt REST API."""

    def __init__(
        self,
        fasta: bytes = UNIPROT_FASTA,
        entry_json: bytes = UNIPROT_JSON,
    ) -> None:
        """Initialize a fake with raw FASTA and complete JSON responses."""
        self.fasta = fasta
        self.entry_json = entry_json
        self.download_calls: list[tuple[str, str]] = []

    def download_entry(self, accession: str, destination: Path) -> None:
        """Record and write only the requested representation."""
        kind, content = {
            ".json": ("entry_json", self.entry_json),
            ".fasta": ("protein_fasta", self.fasta),
        }[destination.suffix]
        self.download_calls.append((kind, accession))
        destination.write_bytes(content)


def _package_files(accession: str = ACCESSION) -> dict[str, bytes]:
    """Return representative files from the requested official include set."""
    accession_root = f"ncbi_dataset/data/{accession}"
    return {
        "README.md": b"NCBI Datasets genome package\n",
        "md5sum.txt": b"official checksums remain untouched\n",
        "ncbi_dataset/data/assembly_data_report.jsonl": (
            json.dumps({"accession": accession}).encode() + b"\n"
        ),
        f"{accession_root}/{FASTA_NAME}": FASTA,
        f"{accession_root}/genomic.gff": b"##gff-version 3\n",
        f"{accession_root}/rna.fna": b">rna\nACGU\n",
        f"{accession_root}/cds_from_genomic.fna": b">cds\nACGT\n",
        f"{accession_root}/protein.faa": b">protein\nMT\n",
        f"{accession_root}/sequence_report.jsonl": b"{}\n",
    }


def _genome_package(
    *,
    accession: str = ACCESSION,
    fasta_path: str | None = None,
    payload: bytes = FASTA,
    catalog_file_count: int = 1,
    include_fasta: bool = True,
    include_gff: bool = True,
    extra_members: dict[str, bytes] | None = None,
) -> bytes:
    """Build a compact, complete NCBI-shaped genome package in memory.

    Returns:
        Complete ZIP package bytes.
    """
    fasta_path = fasta_path or f"{accession}/{FASTA_NAME}"
    catalog_files: list[dict[str, Any]] = [
        {
            "filePath": fasta_path,
            "fileType": "GENOMIC_NUCLEOTIDE_FASTA",
        }
        for _ in range(catalog_file_count)
    ]
    catalog_files.extend(
        [
            {"filePath": f"{accession}/genomic.gff", "fileType": "GFF3"},
            {
                "filePath": f"{accession}/rna.fna",
                "fileType": "RNA_NUCLEOTIDE_FASTA",
            },
            {
                "filePath": f"{accession}/cds_from_genomic.fna",
                "fileType": "CDS_NUCLEOTIDE_FASTA",
            },
            {
                "filePath": f"{accession}/protein.faa",
                "fileType": "PROTEIN_FASTA",
            },
            {
                "filePath": f"{accession}/sequence_report.jsonl",
                "fileType": "SEQUENCE_REPORT",
            },
        ]
    )
    catalog = {
        "apiVersion": "V2",
        "assemblies": [
            {
                "files": [
                    {
                        "filePath": "assembly_data_report.jsonl",
                        "fileType": "DATA_REPORT",
                    }
                ]
            },
            {"accession": accession, "files": catalog_files},
        ],
    }
    package_files = _package_files(accession)
    package_files[f"ncbi_dataset/data/{fasta_path}"] = payload
    if not include_fasta:
        package_files.pop(f"ncbi_dataset/data/{fasta_path}")
    if not include_gff:
        package_files.pop(f"ncbi_dataset/data/{accession}/genomic.gff")
    package_files["ncbi_dataset/data/dataset_catalog.json"] = json.dumps(
        catalog
    ).encode()
    package_files.update(extra_members or {})

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as package:
        for path, content in package_files.items():
            package.writestr(path, content)
    return buffer.getvalue()


def _install_provider(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    downloader: FakeDatasetsCli,
) -> None:
    """Route public path resolution through a test cache and fake CLI."""
    monkeypatch.setattr(artifacts_module.settings, "home", tmp_path)
    monkeypatch.setattr(
        artifacts_module, "download_genome_package", downloader.download_genome_package
    )


def _install_uniprot_provider(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    downloader: FakeUniProtApi,
) -> None:
    """Route UniProt path resolution through a test cache and fake API."""
    monkeypatch.setattr(artifacts_module.settings, "home", tmp_path)
    monkeypatch.setattr(
        artifacts_module, "download_uniprot_entry", downloader.download_entry
    )


@pytest.mark.parametrize(
    "value",
    [
        f"refseq.gcf://{ACCESSION}",
        f"refseq.gcf:{ACCESSION}",
        f"https://identifiers.org/refseq.gcf:{ACCESSION}",
        ACCESSION,
    ],
)
def test_parse_identifier_canonicalizes_each_supported_exact_form(value: str) -> None:
    """Canonicalize URI, Compact, web, and allowlisted bare forms without I/O."""
    identifier = parse_identifier(value)

    assert identifier.namespace == "refseq.gcf"
    assert identifier.accession == ACCESSION
    assert identifier.compact_id == f"refseq.gcf:{ACCESSION}"
    assert identifier.resource_uri == f"refseq.gcf://{ACCESSION}"


@pytest.mark.parametrize(
    "value",
    [
        f"What is the GC content of {ACCESSION}?",
        "P12345",
        "refseq.gcf://GCF_00001030",
        "unknown://value",
        f"{ACCESSION} {ACCESSION}",
    ],
)
def test_parse_identifier_rejects_non_exact_or_non_allowlisted_input(
    value: str,
) -> None:
    """Keep the single-value API separate from permissive prompt scanning."""
    with pytest.raises(IdentifierSyntaxError):
        parse_identifier(value)


def test_parse_identifier_accepts_explicit_uniprot_uri_only() -> None:
    """Recognize the registry URI without adding UniProt to bare-ID inference."""
    identifier = parse_identifier(f"uniprot://{UNIPROT_ACCESSION}")

    assert identifier.namespace == "uniprot"
    assert identifier.accession == UNIPROT_ACCESSION
    assert identifier.compact_id == f"uniprot:{UNIPROT_ACCESSION}"
    with pytest.raises(IdentifierSyntaxError):
        parse_identifier(UNIPROT_ACCESSION)


def test_parse_identifier_uses_generic_resource_for_other_registries() -> None:
    """Canonicalize generic registry IDs without restoring per-registry schemes."""
    identifier = parse_identifier("identifiers://doi:10.1038/s41586-020-2649-2")

    assert identifier.namespace == "doi"
    assert identifier.accession == "10.1038/s41586-020-2649-2"
    assert identifier.compact_id == "doi:10.1038/s41586-020-2649-2"
    assert identifier.resource_uri == ("identifiers://doi:10.1038/s41586-020-2649-2")
    with pytest.raises(IdentifierSyntaxError):
        parse_identifier("doi://10.1038/s41586-020-2649-2")


def test_datasets_cli_uses_the_official_download_command(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Pass the exact requested include set without a shell or REST substitute."""
    commands: list[list[str]] = []

    def fake_run(
        command: list[str],
        *,
        check: bool,
        capture_output: bool,
        text: bool,
        timeout: float,
    ) -> subprocess.CompletedProcess[str]:
        """Record one subprocess invocation and emit its requested ZIP.

        Returns:
            Successful synthetic process result.
        """
        assert check is False
        assert capture_output is True
        assert text is True
        assert timeout == DATASETS_CLI_TIMEOUT_SECONDS
        commands.append(command)
        Path(command[command.index("--filename") + 1]).write_bytes(_genome_package())
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(
        artifacts_module.shutil, "which", lambda _: "/opt/ncbi/datasets"
    )
    monkeypatch.setattr(artifacts_module.subprocess, "run", fake_run)
    destination = tmp_path / "ncbi_dataset.zip"

    download_genome_package(ACCESSION, destination)

    assert DATASETS_INCLUDE == "gff3,rna,cds,protein,genome,seq-report"
    assert commands == [
        [
            "/opt/ncbi/datasets",
            "download",
            "genome",
            "accession",
            ACCESSION,
            "--include",
            DATASETS_INCLUDE,
            "--filename",
            str(destination),
            "--no-progressbar",
        ]
    ]


def test_datasets_cli_returns_official_genome_summary_stdout(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Return the summary command's JSON text without local reshaping."""
    commands: list[list[str]] = []

    def fake_run(
        command: list[str],
        *,
        check: bool,
        capture_output: bool,
        text: bool,
        timeout: float,
    ) -> subprocess.CompletedProcess[str]:
        assert check is False
        assert capture_output is True
        assert text is True
        assert timeout == DATASETS_CLI_TIMEOUT_SECONDS
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, GENOME_SUMMARY, "")

    monkeypatch.setattr(
        artifacts_module.shutil, "which", lambda _: "/opt/ncbi/datasets"
    )
    monkeypatch.setattr(artifacts_module.subprocess, "run", fake_run)

    summary = genome_summary(VERSIONLESS_ACCESSION)

    assert summary == GENOME_SUMMARY
    assert commands == [
        [
            "/opt/ncbi/datasets",
            "summary",
            "genome",
            "accession",
            VERSIONLESS_ACCESSION,
        ]
    ]


@pytest.mark.parametrize(
    ("summary", "error", "message"),
    [
        ("not JSON", ArtifactServiceError, "valid JSON"),
        (
            '{"reports":[],"total_count":0}',
            ArtifactServiceError,
            "exactly one report",
        ),
        (
            '{"reports":[{"accession":"GCF_000006945.1"}],"total_count":1}',
            ArtifactNotFoundError,
            "exact assembly version",
        ),
    ],
)
def test_datasets_cli_rejects_invalid_or_mismatched_genome_summary(
    monkeypatch: pytest.MonkeyPatch,
    summary: str,
    error: type[Exception],
    message: str,
) -> None:
    """Validate the external JSON once at the NCBI CLI boundary."""
    monkeypatch.setattr(
        artifacts_module.shutil, "which", lambda _: "/opt/ncbi/datasets"
    )
    monkeypatch.setattr(
        artifacts_module.subprocess,
        "run",
        lambda command, **_kwargs: subprocess.CompletedProcess(command, 0, summary, ""),
    )

    with pytest.raises(error, match=message):
        genome_summary(ACCESSION)


def test_datasets_cli_missing_or_rejected_is_a_service_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Explain a missing CLI and preserve its diagnostic when it rejects a request."""
    destination = tmp_path / "ncbi_dataset.zip"
    monkeypatch.setattr(artifacts_module.shutil, "which", lambda _: None)
    with pytest.raises(ArtifactServiceError, match=r"install.*datasets"):
        download_genome_package(ACCESSION, destination)

    monkeypatch.setattr(
        artifacts_module.shutil, "which", lambda _: "/opt/ncbi/datasets"
    )
    monkeypatch.setattr(
        artifacts_module.subprocess,
        "run",
        lambda command, **_: subprocess.CompletedProcess(
            command, 1, "", "Error: no genome package"
        ),
    )
    with pytest.raises(ArtifactServiceError, match="no genome package"):
        download_genome_package(ACCESSION, destination)


def test_datasets_cli_timeout_is_a_service_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Bound a hung datasets CLI instead of blocking path resolution forever."""
    destination = tmp_path / "ncbi_dataset.zip"
    monkeypatch.setattr(
        artifacts_module.shutil, "which", lambda _: "/opt/ncbi/datasets"
    )

    def hanging_run(command, **kwargs):
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])

    monkeypatch.setattr(artifacts_module.subprocess, "run", hanging_run)

    with pytest.raises(ArtifactServiceError, match="timed out"):
        download_genome_package(ACCESSION, destination)


def test_uniprot_api_uses_official_entry_urls_and_preserves_responses(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Fetch accession-named JSON and FASTA without parsing or rewriting them."""
    requests: list[tuple[str, str | None, float]] = []

    def fake_urlopen(request: Any, *, timeout: float) -> io.BytesIO:
        """Record the request and return one raw UniProt response.

        Returns:
            In-memory response containing the exact requested bytes.
        """
        requests.append((request.full_url, request.get_header("Accept"), timeout))
        payload = UNIPROT_JSON if request.full_url.endswith(".json") else UNIPROT_FASTA
        return io.BytesIO(payload)

    monkeypatch.setattr(artifacts_module, "urlopen", fake_urlopen)
    fasta_destination = tmp_path / f"{UNIPROT_ACCESSION}.fasta"
    json_destination = tmp_path / f"{UNIPROT_ACCESSION}.json"

    download_uniprot_entry(UNIPROT_ACCESSION, fasta_destination, timeout_seconds=7)
    download_uniprot_entry(UNIPROT_ACCESSION, json_destination, timeout_seconds=7)

    assert requests == [
        (
            f"https://rest.uniprot.org/uniprotkb/{UNIPROT_ACCESSION}.fasta",
            "text/plain",
            7,
        ),
        (
            f"https://rest.uniprot.org/uniprotkb/{UNIPROT_ACCESSION}.json",
            "application/json",
            7,
        ),
    ]
    assert fasta_destination.read_bytes() == UNIPROT_FASTA
    assert json_destination.read_bytes() == UNIPROT_JSON


def test_path_preserves_complete_package_and_skips_cached_cli(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Return the original FASTA path while leaving every package file in place."""
    downloader = FakeDatasetsCli()
    _install_provider(monkeypatch, tmp_path, downloader)

    first = path(ACCESSION)
    second = path(ACCESSION)

    package_root = tmp_path / "artifacts" / "refseq.gcf" / ACCESSION
    extracted_files = {
        member.relative_to(package_root).as_posix()
        for member in package_root.rglob("*")
        if member.is_file()
    }
    expected_files = set(_package_files()) | {"ncbi_dataset/data/dataset_catalog.json"}
    assert extracted_files == expected_files
    assert first.package_root == package_root
    assert first.path == package_root / "ncbi_dataset/data" / ACCESSION / FASTA_NAME
    assert first.path.name == FASTA_NAME
    assert Path(first).read_bytes() == FASTA
    assert os.fspath(first) == str(first.path)
    assert first.requested_identifier.accession == ACCESSION
    assert first.identifier.accession == ACCESSION
    assert first.kind == "genome_fasta"
    assert first.size == len(FASTA)
    assert second == first
    assert downloader.download_calls == [ACCESSION]
    assert not list(package_root.rglob("manifest.json"))
    assert not list(package_root.rglob("genome.*.fna"))


@pytest.mark.parametrize(
    ("kind", "filename", "content"),
    [
        ("genome_fasta", FASTA_NAME, FASTA),
        ("annotation_gff3", "genomic.gff", b"##gff-version 3\n"),
        ("rna_fasta", "rna.fna", b">rna\nACGU\n"),
        ("cds_fasta", "cds_from_genomic.fna", b">cds\nACGT\n"),
        ("protein_fasta", "protein.faa", b">protein\nMT\n"),
    ],
)
def test_path_selects_each_refseq_catalog_kind(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    kind: str,
    filename: str,
    content: bytes,
) -> None:
    """Select the original catalog member matching each registered kind."""
    downloader = FakeDatasetsCli()
    _install_provider(monkeypatch, tmp_path, downloader)

    artifact = path(ACCESSION, artifact=kind)

    package_root = tmp_path / "artifacts" / "refseq.gcf" / ACCESSION
    assert artifact.path == package_root / "ncbi_dataset/data" / ACCESSION / filename
    assert artifact.path.read_bytes() == content
    assert artifact.kind == kind
    assert artifact.size == len(content)
    assert downloader.download_calls == [ACCESSION]


def test_one_cached_package_serves_every_refseq_kind(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Serve all RefSeq kinds from one downloaded complete package."""
    downloader = FakeDatasetsCli()
    _install_provider(monkeypatch, tmp_path, downloader)

    genome = path(ACCESSION)
    annotation = path(ACCESSION, artifact="annotation_gff3")
    protein = path(ACCESSION, artifact="protein_fasta")

    assert genome.kind == "genome_fasta"
    assert annotation.path.name == "genomic.gff"
    assert protein.path.name == "protein.faa"
    assert genome.package_root == annotation.package_root == protein.package_root
    assert downloader.download_calls == [ACCESSION]


def test_path_rejects_package_missing_the_requested_kind(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Report a stable package error when the catalog member is absent."""
    downloader = FakeDatasetsCli(_genome_package(include_gff=False))
    _install_provider(monkeypatch, tmp_path, downloader)

    with pytest.raises(ArtifactPackageError, match="missing"):
        path(ACCESSION, artifact="annotation_gff3")

    assert not (tmp_path / "artifacts" / "refseq.gcf" / ACCESSION).exists()


def test_versionless_path_uses_catalog_canonical_accession(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Let the official package identify the version selected by NCBI."""
    downloader = FakeDatasetsCli()
    _install_provider(monkeypatch, tmp_path, downloader)

    first = path(VERSIONLESS_ACCESSION)
    second = path(VERSIONLESS_ACCESSION)

    assert first.requested_identifier.accession == VERSIONLESS_ACCESSION
    assert first.identifier.accession == ACCESSION
    assert second.identifier.accession == ACCESSION
    assert downloader.download_calls == [VERSIONLESS_ACCESSION]


def test_versioned_path_never_silently_upgrades(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Reject a package for a different version than the explicit request."""
    downloader = FakeDatasetsCli(_genome_package(accession="GCF_000006945.3"))
    _install_provider(monkeypatch, tmp_path, downloader)

    with pytest.raises(ArtifactNotFoundError, match="exact assembly version"):
        path(ACCESSION)


@pytest.mark.parametrize(
    ("package", "message"),
    [
        (
            _genome_package(extra_members={"../escape.txt": b"bad"}),
            "unsafe",
        ),
        (_genome_package(extra_members={"nested//file.txt": b"bad"}), "unsafe"),
        (_genome_package(extra_members={"C:/escape.txt": b"bad"}), "unsafe"),
        (_genome_package(extra_members={"empty//": b""}), "unsafe"),
        (_genome_package(catalog_file_count=2), "exactly one"),
        (_genome_package(include_fasta=False), "missing"),
        (b"not a zip", "valid ZIP"),
    ],
)
def test_path_rejects_unsafe_or_invalid_complete_packages(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    package: bytes,
    message: str,
) -> None:
    """Reject invalid packages without publishing partial extracted content."""
    downloader = FakeDatasetsCli(package)
    _install_provider(monkeypatch, tmp_path, downloader)

    with pytest.raises(ArtifactPackageError, match=message):
        path(ACCESSION)

    cache_path = tmp_path / "artifacts" / "refseq.gcf" / ACCESSION
    assert not cache_path.exists()
    assert not (tmp_path / "artifacts" / "escape.txt").exists()


def test_failed_package_does_not_poison_later_cache_fill(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Leave no cache hit behind after an invalid package download."""
    downloader = FakeDatasetsCli(b"not a zip")
    _install_provider(monkeypatch, tmp_path, downloader)

    with pytest.raises(ArtifactPackageError):
        path(ACCESSION)
    downloader.package = _genome_package()

    artifact = path(ACCESSION)

    assert Path(artifact).read_bytes() == FASTA
    assert downloader.download_calls == [ACCESSION, ACCESSION]


def test_uniprot_json_and_fasta_cache_independently_on_demand(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Fetch only the requested raw representation and reuse each cached file."""
    downloader = FakeUniProtApi()
    _install_uniprot_provider(monkeypatch, tmp_path, downloader)

    metadata = path(
        f"uniprot:{UNIPROT_ACCESSION}",
        artifact="entry_json",
    )
    cached_metadata = path(
        f"uniprot:{UNIPROT_ACCESSION}",
        artifact="entry_json",
    )
    package_root = tmp_path / "artifacts" / "uniprot" / UNIPROT_ACCESSION
    json_path = package_root / f"{UNIPROT_ACCESSION}.json"

    assert metadata.path == json_path
    assert metadata.path.read_bytes() == UNIPROT_JSON
    assert metadata.kind == "entry_json"
    assert metadata.size == len(UNIPROT_JSON)
    assert cached_metadata == metadata
    assert downloader.download_calls == [("entry_json", UNIPROT_ACCESSION)]
    assert [member.name for member in package_root.iterdir()] == [
        f"{UNIPROT_ACCESSION}.json"
    ]

    first = path(f"uniprot://{UNIPROT_ACCESSION}")
    second = path(f"uniprot:{UNIPROT_ACCESSION}")

    assert first.path == package_root / f"{UNIPROT_ACCESSION}.fasta"
    assert first.package_root == package_root
    assert first.path.read_bytes() == UNIPROT_FASTA
    assert first.requested_identifier.accession == UNIPROT_ACCESSION
    assert first.identifier.accession == UNIPROT_ACCESSION
    assert first.kind == "protein_fasta"
    assert first.size == len(UNIPROT_FASTA)
    assert second == first
    assert downloader.download_calls == [
        ("entry_json", UNIPROT_ACCESSION),
        ("protein_fasta", UNIPROT_ACCESSION),
    ]
    assert sorted(member.name for member in package_root.iterdir()) == [
        f"{UNIPROT_ACCESSION}.fasta",
        f"{UNIPROT_ACCESSION}.json",
    ]
    entry = json.loads(json_path.read_bytes())
    pdb_ids = [
        reference["id"]
        for reference in entry["uniProtKBCrossReferences"]
        if reference["database"] == "PDB"
    ]
    assert pdb_ids == ["1EMA", "1GFL"]


def test_uniprot_fasta_cache_does_not_require_json(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Reuse FASTA independently and fetch JSON only when explicitly requested."""
    downloader = FakeUniProtApi()
    _install_uniprot_provider(monkeypatch, tmp_path, downloader)
    package_root = tmp_path / "artifacts" / "uniprot" / UNIPROT_ACCESSION
    package_root.mkdir(parents=True)
    (package_root / f"{UNIPROT_ACCESSION}.fasta").write_bytes(UNIPROT_FASTA)

    fasta = path(f"uniprot://{UNIPROT_ACCESSION}")
    metadata = path(
        f"uniprot://{UNIPROT_ACCESSION}",
        artifact="entry_json",
    )

    assert fasta.path.read_bytes() == UNIPROT_FASTA
    assert metadata.path.read_bytes() == UNIPROT_JSON
    assert downloader.download_calls == [("entry_json", UNIPROT_ACCESSION)]


def test_invalid_uniprot_response_does_not_poison_later_cache_fill(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Publish only a FASTA whose UniProt header matches the requested entry."""
    downloader = FakeUniProtApi(b">sp|P12345|WRONG_ENTRY Wrong protein\nMA\n")
    _install_uniprot_provider(monkeypatch, tmp_path, downloader)

    with pytest.raises(ArtifactNotFoundError, match="exact accession"):
        path(f"uniprot://{UNIPROT_ACCESSION}")
    downloader.fasta = UNIPROT_FASTA

    artifact = path(f"uniprot://{UNIPROT_ACCESSION}")

    assert artifact.path.read_bytes() == UNIPROT_FASTA
    assert downloader.download_calls == [
        ("protein_fasta", UNIPROT_ACCESSION),
        ("protein_fasta", UNIPROT_ACCESSION),
    ]


def test_invalid_uniprot_json_does_not_poison_later_cache_fill(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Publish only JSON whose primary accession matches the requested entry."""
    downloader = FakeUniProtApi(entry_json=b'{"primaryAccession":"P12345"}\n')
    _install_uniprot_provider(monkeypatch, tmp_path, downloader)

    with pytest.raises(ArtifactNotFoundError, match="exact accession"):
        path(
            f"uniprot://{UNIPROT_ACCESSION}",
            artifact="entry_json",
        )
    downloader.entry_json = UNIPROT_JSON

    artifact = path(
        f"uniprot://{UNIPROT_ACCESSION}",
        artifact="entry_json",
    )

    assert artifact.path.read_bytes() == UNIPROT_JSON
    assert downloader.download_calls == [
        ("entry_json", UNIPROT_ACCESSION),
        ("entry_json", UNIPROT_ACCESSION),
    ]


@pytest.mark.parametrize("mode", ["rb", "rt"])
def test_open_works_as_thin_standard_file_adapter(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    mode: str,
) -> None:
    """Open the original cached FASTA without an analysis wrapper."""
    downloader = FakeDatasetsCli()
    _install_provider(monkeypatch, tmp_path, downloader)

    with open(ACCESSION, mode=mode) as fasta_file:
        assert fasta_file.read() == (FASTA if mode == "rb" else FASTA.decode())


def test_path_rejects_unsupported_namespace_and_artifact_kind() -> None:
    """Dispatch only explicitly registered namespace and artifact pairs."""
    with pytest.raises(UnsupportedArtifactError, match="uniprot"):
        path("uniprot:P12345", artifact="genome_fasta")
    with pytest.raises(UnsupportedArtifactError, match="entry_json"):
        path(ACCESSION, artifact="entry_json")
