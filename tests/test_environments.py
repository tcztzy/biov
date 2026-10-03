"""Fresh-device Pixi setup, integrity checks, and locked manifest installation."""

import errno
import hashlib
import io
import json
import shutil
import subprocess
import tarfile
import tomllib
import zipfile
from pathlib import Path

import pytest
from typer.testing import CliRunner

from biov import environments
from biov.cli import app
from biov.config import Settings


@pytest.fixture(autouse=True)
def managed_pixi_only(monkeypatch):
    """Resolve the managed copy unless a test provides its own Pixi candidate."""
    monkeypatch.delenv(environments.PIXI_BIN, raising=False)
    monkeypatch.setattr(environments.shutil, "which", lambda name: None)


def manifest_environments(manifest: Path) -> list[str]:
    """Return the environment names one manifest declares.

    Returns:
        Declared names, in declaration order.
    """
    with manifest.open("rb") as source:
        document = tomllib.load(source)
    records = document.get("tool", {}).get("pixi", {}).get("environments", {})
    return list(records) if isinstance(records, dict) else []


def environment_prefix(manifest: Path, name: str) -> Path:
    """Return the prefix Pixi installs one environment into.

    Returns:
        Prefix path under the manifest's own project directory.
    """
    return manifest.parent / ".pixi" / "envs" / name


def pixi_info_response(argv: list[str]) -> str:
    """Return the ``pixi info --json`` document a faked manager reports.

    Returns:
        JSON text naming each declared environment's prefix.
    """
    manifest = Path(argv[argv.index("--manifest-path") + 1])
    return json.dumps(
        {
            "environments_info": [
                {"name": name, "prefix": str(environment_prefix(manifest, name))}
                for name in manifest_environments(manifest)
            ]
        }
    )


def fake_pixi_run(argv: list[str]) -> subprocess.CompletedProcess[str]:
    """Answer one faked Pixi invocation and create what an install would create.

    Returns:
        Successful synthetic process result.
    """
    if argv[1] == "info":
        return subprocess.CompletedProcess(argv, 0, pixi_info_response(argv))
    if argv[1] == "install":
        manifest = Path(argv[argv.index("--manifest-path") + 1])
        names = (
            manifest_environments(manifest)
            if "--all" in argv
            else [argv[argv.index("--environment") + 1]]
        )
        for name in names:
            environment_prefix(manifest, name).mkdir(parents=True, exist_ok=True)
    return subprocess.CompletedProcess(argv, 0, f"pixi {environments.PIXI_VERSION}\n")


def register_release(monkeypatch, directory: Path) -> Path:
    """Register a checksum-verified release archive for a faked Linux host.

    Returns:
        Path of the release archive the faked download would serve.
    """
    archive = directory / "pixi.tar.gz"
    with tarfile.open(archive, "w:gz") as package:
        entry = tarfile.TarInfo("pixi")
        entry.size = len(b"executable")
        entry.mode = 0o755
        package.addfile(entry, io.BytesIO(b"executable"))
    monkeypatch.setattr(environments.platform, "system", lambda: "Linux")
    monkeypatch.setattr(environments.platform, "machine", lambda: "x86_64")
    monkeypatch.setattr(
        environments,
        "_RELEASES",
        {
            ("Linux", "x86_64"): (
                archive.name,
                hashlib.sha256(archive.read_bytes()).hexdigest(),
            )
        },
    )
    return archive


def report_pinned_version(monkeypatch, calls: list[list[str]] | None = None):
    """Answer every Pixi invocation with the pinned version and record argv."""

    def run(argv, **kwargs):
        if calls is not None:
            calls.append(argv)
        return fake_pixi_run(argv)

    monkeypatch.setattr(environments.subprocess, "run", run)


@pytest.mark.parametrize("suffix", ["tar.gz", "zip"])
def test_setup_installs_verified_manager_and_locked_environment(
    monkeypatch, tmp_path, suffix
):
    """Reject corrupt bundles, retain notices, reuse installs, and require the lock."""
    archive = tmp_path / f"pixi.{suffix}"
    if suffix == "zip":
        with zipfile.ZipFile(archive, "w") as package:
            package.writestr("pixi", b"executable")
    else:
        with tarfile.open(archive, "w:gz") as package:
            member = tarfile.TarInfo("pixi")
            member.size = len(b"executable")
            member.mode = 0o755
            package.addfile(member, io.BytesIO(b"executable"))
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    monkeypatch.setattr(environments.platform, "system", lambda: "Linux")
    monkeypatch.setattr(environments.platform, "machine", lambda: "x86_64")
    monkeypatch.setattr(
        environments, "_RELEASES", {("Linux", "x86_64"): (archive.name, digest)}
    )
    config = Settings.model_validate(
        {
            "environment_root": tmp_path / "environments",
        }
    )
    monkeypatch.setattr(environments, "settings", config)
    executable = environments.pixi_path()
    calls = []

    def run(argv, **kwargs):
        assert kwargs == (
            {"check": True, "capture_output": True, "text": True}
            if kwargs.get("capture_output")
            else {"check": True}
        )
        calls.append(argv)
        return subprocess.CompletedProcess(
            argv, 0, f"pixi {environments.PIXI_VERSION}\n"
        )

    monkeypatch.setattr(environments.subprocess, "run", run)
    bad = tmp_path / "bad.archive"
    bad.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="SHA-256"):
        environments.setup_environment(archive=bad)
    assert not executable.exists()
    assert not calls
    assert environments.setup_environment(archive=archive) == executable
    assert "BSD 3-Clause" in (executable.parent / "LICENSE").read_text()
    assert not environments.environment_manifest().exists()
    archive.unlink()
    calls.clear()
    manifest = environments.setup_environment("default")
    assert manifest.is_file()
    assert manifest.with_name("pixi.lock").is_file()
    assert calls[1] == [
        str(executable),
        "install",
        "--no-config",
        "--manifest-path",
        str(manifest),
        "--locked",
        "--environment",
        "default",
    ]
    original = manifest.read_bytes()
    assert environments.setup_environment("default") == manifest
    assert manifest.read_bytes() == original
    # A caller-owned project is used in place, never replaced by the bundled one.
    project = tmp_path / "pyproject.toml"
    project.write_text(
        'description = "custom project"\n'
        "[tool.pixi.environments]\n"
        "default = { features = [], no-default-feature = true }\n"
    )
    config.environment_manifest = project
    assert environments.setup_environment("default") == project
    assert "custom project" in project.read_text()

    def failed_install(argv, **kwargs):
        raise subprocess.CalledProcessError(1, argv)

    monkeypatch.setattr(environments.subprocess, "run", failed_install)
    with pytest.raises(subprocess.CalledProcessError):
        environments.setup_environment("default")


def test_setup_download_failure_and_remote_rejection(monkeypatch, tmp_path):
    """Fail without installing another runtime or hiding failed downloads."""
    config = Settings(environment_root=tmp_path / "environments")
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(environments.platform, "system", lambda: "Linux")
    monkeypatch.setattr(environments.platform, "machine", lambda: "x86_64")
    urls = []

    def fail(url, **kwargs):
        urls.append(url)
        raise TimeoutError("download stopped")

    monkeypatch.setattr(environments, "urlopen", fail)
    with pytest.raises(TimeoutError):
        environments.setup_environment()
    assert urls == [
        f"https://github.com/prefix-dev/pixi/releases/download/v{environments.PIXI_VERSION}/pixi-x86_64-unknown-linux-musl.tar.gz"
    ]
    assert not environments.pixi_path().exists()
    assert not list(config.environment_root.glob(".pixi-*"))
    config.execution_host = "compute"
    with pytest.raises(ValueError, match="execution host"):
        environments.setup_environment()
    assert len(urls) == 1


def test_declared_environments_resolve_and_bulk_setup_prepares(monkeypatch, tmp_path):
    """Install declared environments and run their preparation tasks."""
    monkeypatch.setattr(environments.platform, "system", lambda: "Linux")
    monkeypatch.setattr(environments.platform, "machine", lambda: "x86_64")
    config = Settings(environment_root=tmp_path)
    monkeypatch.setattr(environments, "settings", config)
    declared = environments.declared_environments()
    assert {"default", "python", "samtools", "uce"} <= set(declared)
    assert environments.resolve_environment("python") == "python"
    assert environments.resolve_environment("samtools") == "samtools"
    assert environments.resolve_environment("uce") == "uce"
    executable = environments.pixi_path()
    executable.parent.mkdir(parents=True)
    executable.touch()
    calls = []
    monkeypatch.setattr(
        environments.subprocess,
        "run",
        lambda argv, **kw: (
            calls.append(argv)
            or subprocess.CompletedProcess(
                argv, 0, f"pixi {environments.PIXI_VERSION}\n"
            )
        ),
    )
    manifest = environments.setup_environment(all_environments=True)
    assert calls[1] == [
        str(executable),
        "install",
        "--no-config",
        "--manifest-path",
        str(manifest),
        "--locked",
        "--all",
    ]
    # Only environments that declare a prepare task are prepared, in the order
    # the manifest declares them.
    assert calls[2] == [
        str(executable),
        "run",
        "--no-config",
        "--manifest-path",
        str(manifest),
        "--as-is",
        "--environment",
        "diffdock",
        "prepare",
    ]
    assert [
        argv[argv.index("--environment") + 1] for argv in calls if argv[1] == "run"
    ] == ["diffdock", "autosite", "r", "uce"]
    with pytest.raises(ValueError, match="one environment"):
        environments.setup_environment("samtools", all_environments=True)
    with pytest.raises(ValueError, match="BIOV_ENVIRONMENT_MANIFEST"):
        environments.setup_environment(update_lock=True)
    config.environment_manifest = manifest
    calls.clear()
    environments.setup_environment(update_lock=True)
    assert calls[1] == [
        str(executable),
        "lock",
        "--no-config",
        "--manifest-path",
        str(manifest),
    ]


def test_runtime_manifest_declares_environments_and_prepare_tasks(
    monkeypatch, tmp_path
):
    """Ship environments, features and prepare tasks, never command mappings."""
    import tomllib

    config = Settings(environment_root=tmp_path)
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(environments.platform, "system", lambda: "Linux")
    monkeypatch.setattr(environments.platform, "machine", lambda: "x86_64")
    document = tomllib.loads(environments.manifest_source().read_text())["tool"]
    assert set(document) == {"pixi"}
    pixi = document["pixi"]
    declared = environments.declared_environments()
    assert set(pixi["feature"]) == set(declared) - {"default"}
    for name, features in declared.items():
        assert pixi["environments"][name] == {
            "features": list(features),
            "no-default-feature": True,
        }
        assert set(features) <= set(pixi["feature"])
    # Setup prepares exactly the environments whose features declare the task.
    assert {
        name
        for name, record in pixi["feature"].items()
        if environments.PREPARE_TASK in record.get("tasks", {})
    } == {"autosite", "diffdock", "r", "uce"}
    assert environments.prepare_environments() == ("diffdock", "autosite", "r", "uce")
    # Every prepared environment resolves to a real task, and the ones that ship
    # no preparation resolve to none at all.
    for name in environments.prepare_environments():
        assert environments.prepare_command(name)
    assert environments.prepare_command("samtools") is None
    executable = environments.pixi_path()
    executable.parent.mkdir(parents=True)
    executable.touch()
    calls = []

    def run(argv, **kwargs):
        calls.append(argv)
        return subprocess.CompletedProcess(
            argv, 0, f"pixi {environments.PIXI_VERSION}\n"
        )

    monkeypatch.setattr(environments.subprocess, "run", run)
    installed = environments.setup_environment("uce")
    assert calls[1] == [
        str(executable),
        "install",
        "--no-config",
        "--manifest-path",
        str(installed),
        "--locked",
        "--environment",
        "uce",
    ]
    assert calls[2] == [
        str(executable),
        "run",
        "--no-config",
        "--manifest-path",
        str(installed),
        "--as-is",
        "--environment",
        "uce",
        "prepare",
    ]
    calls.clear()
    environments.setup_environment("samtools")
    assert [argv[1] for argv in calls] == ["--version", "install"]
    monkeypatch.setattr(
        environments.subprocess,
        "run",
        lambda *a, **kw: subprocess.CompletedProcess(a, 0, "pixi 9.9.9"),
    )
    with pytest.raises(ValueError, match="Expected pixi"):
        environments.setup_environment()


def test_workspace_publication_race(monkeypatch, tmp_path):
    """Reuse a complete workspace published by another setup process."""
    monkeypatch.setattr(environments, "settings", Settings(environment_root=tmp_path))

    def competing_publish(source, target):
        shutil.copytree(source, target)
        raise OSError(errno.ENOTEMPTY, "another process published first")

    monkeypatch.setattr(environments.os, "rename", competing_publish)
    manifest = environments.environment_manifest(prepare=True)
    assert manifest.read_bytes() == environments.manifest_source().read_bytes()
    assert manifest.with_name("pixi.lock").is_file()


def test_pixi_tar_rejects_external_links(monkeypatch, tmp_path):
    """Apply the data filter even on Python versions whose default is permissive."""
    archive = tmp_path / "hostile.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        member = tarfile.TarInfo("pixi")
        member.type = tarfile.SYMTYPE
        member.linkname = str(tmp_path / "outside")
        bundle.addfile(member)
    monkeypatch.setattr(
        environments, "settings", Settings(environment_root=tmp_path / "manager")
    )
    monkeypatch.setattr(environments.platform, "system", lambda: "Linux")
    monkeypatch.setattr(environments.platform, "machine", lambda: "x86_64")
    monkeypatch.setattr(
        environments,
        "_RELEASES",
        {
            ("Linux", "x86_64"): (
                archive.name,
                hashlib.sha256(archive.read_bytes()).hexdigest(),
            )
        },
    )
    extractall = tarfile.TarFile.extractall

    def filtered(self, *args, **kwargs):
        assert kwargs["filter"] == "data"
        return extractall(self, *args, **kwargs)

    monkeypatch.setattr(tarfile.TarFile, "extractall", filtered)
    with pytest.raises(tarfile.AbsoluteLinkError):
        environments.setup_environment(archive=archive)
    assert not environments.pixi_path().exists()


def test_matching_system_pixi_is_reused_without_installing(monkeypatch, tmp_path):
    """Use a PATH Pixi of the pinned version instead of downloading a second copy."""
    config = Settings(environment_root=tmp_path / "environments")
    monkeypatch.setattr(environments, "settings", config)
    system = tmp_path / "bin" / "pixi"
    system.parent.mkdir()
    system.write_bytes(b"executable")
    monkeypatch.setattr(environments.shutil, "which", lambda name: str(system))
    calls = []
    report_pinned_version(monkeypatch, calls)
    assert environments.pixi_path() == system
    assert environments.pixi_command() == (str(system),)
    assert environments.setup_environment() == system
    assert calls and all(argv == [str(system), "--version"] for argv in calls)
    assert not environments.managed_pixi_path().exists()
    assert not config.environment_root.exists()


def test_mismatched_system_pixi_warns_and_installs_pinned_copy(monkeypatch, tmp_path):
    """Report a different Pixi on PATH, then install the pinned release."""
    archive = register_release(monkeypatch, tmp_path)
    config = Settings(environment_root=tmp_path / "environments")
    monkeypatch.setattr(environments, "settings", config)
    system = tmp_path / "bin" / "pixi"
    system.parent.mkdir()
    system.write_bytes(b"another pixi")
    monkeypatch.setattr(environments.shutil, "which", lambda name: str(system))
    calls = []

    def run(argv, **kwargs):
        calls.append(argv)
        reported = (
            "pixi 0.0.1\n"
            if argv[0] == str(system)
            else f"pixi {environments.PIXI_VERSION}\n"
        )
        return subprocess.CompletedProcess(argv, 0, reported)

    monkeypatch.setattr(environments.subprocess, "run", run)
    with pytest.warns(RuntimeWarning, match="Ignoring pixi"):
        installed = environments.setup_environment(archive=archive)
    assert installed == environments.managed_pixi_path()
    assert installed.is_file()
    assert (installed.parent / "LICENSE").is_file()
    assert system.read_bytes() == b"another pixi"
    assert calls[0] == [str(system), "--version"]


def test_explicit_pixi_override_wins_and_must_be_usable(monkeypatch, tmp_path):
    """An operator override outranks PATH, and a broken override fails loudly."""
    config = Settings(environment_root=tmp_path / "environments")
    monkeypatch.setattr(environments, "settings", config)
    override = tmp_path / "custom" / "pixi"
    system = tmp_path / "bin" / "pixi"
    monkeypatch.setenv(environments.PIXI_BIN, str(override))
    monkeypatch.setattr(
        environments.shutil,
        "which",
        lambda name: str(system) if name == "pixi" else None,
    )
    assert environments.pixi_path() == override
    assert environments.pixi_command() == (str(override),)
    with pytest.raises(ValueError, match="BIOV_PIXI_BIN"):
        environments.setup_environment()
    override.parent.mkdir()
    override.write_bytes(b"executable")
    calls = []
    report_pinned_version(monkeypatch, calls)
    assert environments.setup_environment() == override
    assert calls and all(argv == [str(override), "--version"] for argv in calls)
    assert not config.environment_root.exists()
    # A bare command name is resolved on PATH like any other executable.
    monkeypatch.setenv(environments.PIXI_BIN, "pixi")
    assert environments.pixi_path() == system


def test_override_with_the_wrong_version_is_rejected(monkeypatch, tmp_path):
    """Never silently accept an override that is not the pinned Pixi."""
    monkeypatch.setattr(
        environments, "settings", Settings(environment_root=tmp_path / "environments")
    )
    override = tmp_path / "custom" / "pixi"
    override.parent.mkdir()
    override.write_bytes(b"executable")
    monkeypatch.setenv(environments.PIXI_BIN, str(override))
    monkeypatch.setattr(
        environments.subprocess,
        "run",
        lambda *a, **kw: subprocess.CompletedProcess(a, 0, "pixi 9.9.9\n"),
    )
    with pytest.raises(ValueError, match="Expected pixi") as failure:
        environments.setup_environment()
    assert str(override) in str(failure.value)


def test_setup_repairs_a_stale_managed_directory(monkeypatch, tmp_path):
    """Clear a managed directory whose binary is gone, then publish and reuse it."""
    archive = register_release(monkeypatch, tmp_path)
    config = Settings(environment_root=tmp_path / "environments")
    monkeypatch.setattr(environments, "settings", config)
    report_pinned_version(monkeypatch)
    executable = environments.managed_pixi_path()
    executable.parent.mkdir(parents=True)
    (executable.parent / "partial.tmp").write_bytes(b"leftover")
    assert not executable.exists()
    assert environments.setup_environment(archive=archive) == executable
    assert executable.is_file()
    assert not (executable.parent / "partial.tmp").exists()
    archive.unlink()
    calls = []
    report_pinned_version(monkeypatch, calls)
    assert environments.setup_environment(archive=archive) == executable
    assert calls == [[str(executable), "--version"]]


def test_manager_publication_race(monkeypatch, tmp_path):
    """Keep a manager published by another setup process."""
    archive = register_release(monkeypatch, tmp_path)
    monkeypatch.setattr(
        environments, "settings", Settings(environment_root=tmp_path / "environments")
    )
    report_pinned_version(monkeypatch)
    executable = environments.managed_pixi_path()

    def competing_publish(source, target):
        shutil.copytree(source, target)
        raise OSError(errno.ENOTEMPTY, "another process published first")

    monkeypatch.setattr(environments.os, "rename", competing_publish)
    assert environments.setup_environment(archive=archive) == executable
    assert executable.is_file()


def test_unknown_environment_names_the_declared_choices(monkeypatch, tmp_path):
    """Fail once, naming the valid environments, when setup cannot resolve itself."""
    config = Settings(environment_root=tmp_path)
    monkeypatch.setattr(environments, "settings", config)
    with pytest.raises(ValueError) as failure:
        environments.setup_environment("nope")
    message = str(failure.value)
    assert "nope" in message
    assert str(environments.manifest_source()) in message
    assert "samtools" in message and "uce" in message
    # The operator sees that same single error through the CLI, not a traceback.
    result = CliRunner().invoke(app, ["setup", "nope"])
    assert result.exit_code == 2, result.output
    assert "nope" in result.output and "samtools" in result.output
    assert "Traceback" not in result.output


def test_native_and_table_environment_features(monkeypatch, tmp_path):
    """Prepare inherited default tasks and feature overrides, honoring opt-outs."""
    manifest = tmp_path / "pyproject.toml"
    manifest.write_text(
        "[tool.pixi.environments]\n"
        'science = ["science"]\n'
        "empty = []\n"
        'analysis = { features = ["science", "analysis"], no-default-feature = true }\n'
        "skipped = { features = [], no-default-feature = true }\n"
        "inherited = { features = [] }\n"
        "[tool.pixi.environments.table]\n"
        'features = ["analysis", "science"]\n'
        "no-default-feature = true\n"
        "[tool.pixi.tasks]\n"
        'prepare = "echo default"\n'
        "[tool.pixi.feature.science.tasks]\n"
        'prepare = "echo ready"\n'
    )
    monkeypatch.setattr(
        environments,
        "settings",
        Settings(
            environment_root=tmp_path / "environments", environment_manifest=manifest
        ),
    )
    assert environments.declared_environments() == {
        "science": ("science",),
        "empty": (),
        "analysis": ("science", "analysis"),
        "skipped": (),
        "inherited": (),
        "table": ("analysis", "science"),
    }
    prepared = ("science", "empty", "analysis", "inherited", "table")
    assert environments.prepare_environments() == prepared
    assert environments.resolve_environment("science") == "science"
    executable = environments.managed_pixi_path()
    executable.parent.mkdir(parents=True)
    executable.touch()
    calls = []
    report_pinned_version(monkeypatch, calls)
    assert environments.setup_environment(all_environments=True) == manifest
    # Resolve the selected task through Pixi, so feature tasks override defaults.
    assert [argv[-2:] for argv in calls if argv[1] == "run"] == [
        [name, "prepare"] for name in prepared
    ]


@pytest.mark.parametrize(
    "declaration",
    [
        'science = "python"',
        "science = 42",
        'science = ["science", 1]',
        'science = { features = "science" }',
        'science = { features = ["science", 1] }',
        "[tool.pixi.environments.science]\nfeatures = 1",
    ],
)
def test_malformed_declared_environment_is_reported(monkeypatch, tmp_path, declaration):
    """Report malformed shorthand and table feature declarations by name."""
    manifest = tmp_path / "pyproject.toml"
    manifest.write_text(f"[tool.pixi.environments]\n{declaration}\n")
    monkeypatch.setattr(
        environments, "settings", Settings(environment_manifest=manifest)
    )
    with pytest.raises(ValueError, match=r"tool\.pixi\.environments\.science"):
        environments.declared_environments()


def test_prepare_command_and_record_invalidate_themselves(monkeypatch, tmp_path):
    """Track preparation per environment, prepare task, lock and prefix."""
    manifest = tmp_path / "pyproject.toml"

    def declare(task: str) -> None:
        manifest.write_text(
            "[tool.pixi.environments]\n"
            'science = { features = ["science"], no-default-feature = true }\n'
            "skipped = { features = [], no-default-feature = true }\n"
            "[tool.pixi.feature.science.tasks]\n"
            f'prepare = "{task}"\n'
        )

    declare("echo ready")
    manifest.with_name("pixi.lock").write_text("lock = 1\n")
    monkeypatch.setattr(
        environments,
        "settings",
        Settings(
            environment_root=tmp_path / "environments", environment_manifest=manifest
        ),
    )
    assert environments.prepare_environments() == ("science",)
    assert environments.prepare_command("science") == json.dumps(["echo ready"])
    # Opting out of the default feature, and naming no environment at all, both
    # resolve to no preparation instead of the workspace default task.
    assert environments.prepare_command("skipped") is None
    assert environments.prepare_command("unlisted") is None
    prefix = environment_prefix(manifest, "science")
    prefix.mkdir(parents=True)
    marker = environments.prepared_marker(prefix)
    assert marker == prefix / environments.PREPARED_RECORD
    assert not environments.preparation_recorded(manifest, "science", prefix)
    environments.record_preparation(manifest, "science", prefix)
    assert marker.is_file()
    assert environments.preparation_recorded(manifest, "science", prefix)
    # Editing the prepare task invalidates the record.
    declare("echo changed")
    assert not environments.preparation_recorded(manifest, "science", prefix)
    # So does a lock the environment was not prepared from.
    declare("echo ready")
    assert environments.preparation_recorded(manifest, "science", prefix)
    manifest.with_name("pixi.lock").write_text("lock = 2\n")
    assert not environments.preparation_recorded(manifest, "science", prefix)
    # Only one truth is left: the record lives in the environment it prepared,
    # so Pixi cleaning or reinstalling that prefix invalidates it too.
    manifest.with_name("pixi.lock").write_text("lock = 1\n")
    declare("echo ready")
    assert environments.preparation_recorded(manifest, "science", prefix)
    shutil.rmtree(prefix)
    assert not environments.preparation_recorded(manifest, "science", prefix)
    assert not marker.exists()
    assert not list(manifest.parent.glob(".biov*"))
    # An environment Pixi cannot report records nothing instead of skipping.
    environments.record_preparation(manifest, "science", None)
    assert not environments.preparation_recorded(manifest, "science", None)


def test_setup_records_preparation_for_on_demand_exec(monkeypatch, tmp_path):
    """Explicit setup still prepares, and records it for the next exec."""
    manifest = tmp_path / "pyproject.toml"
    manifest.write_text(
        "[tool.pixi.environments]\n"
        'uce = { features = ["uce"], no-default-feature = true }\n'
        "[tool.pixi.feature.uce.tasks]\n"
        'prepare = "git fetch uce"\n'
    )
    manifest.with_name("pixi.lock").write_text("lock = 1\n")
    monkeypatch.setattr(
        environments,
        "settings",
        Settings(
            environment_root=tmp_path / "environments", environment_manifest=manifest
        ),
    )
    executable = environments.managed_pixi_path()
    executable.parent.mkdir(parents=True)
    executable.touch()
    calls = []
    report_pinned_version(monkeypatch, calls)
    assert environments.setup_environment("uce") == manifest
    assert [argv[1] for argv in calls[-2:]] == ["run", "info"]
    assert calls[-2][-1] == "prepare"
    prefixes = environments.environment_prefixes(str(executable), manifest)
    assert environments.preparation_recorded(manifest, "uce", prefixes["uce"])
    # Setup prepares unconditionally on every explicit call; an exec after it
    # installs and runs without repeating minutes of preparation.
    calls.clear()
    environments.setup_environment("uce")
    assert [argv[1] for argv in calls if argv[1] != "info"] == [
        "--version",
        "install",
        "run",
    ]
    calls.clear()
    environments.provision_environment("uce")
    assert [argv[1] for argv in calls] == ["install", "info"]
    # Pixi cleaning the environment deletes the prefix and its record together,
    # so the next exec prepares the reinstalled environment again.
    shutil.rmtree(environment_prefix(manifest, "uce"))
    assert not environments.preparation_recorded(
        manifest, "uce", environment_prefix(manifest, "uce")
    )
    calls.clear()
    environments.provision_environment("uce")
    assert [argv[1] for argv in calls] == ["install", "info", "run"]
    assert calls[-1][-1] == "prepare"


def test_platform_scoped_prepare_tasks_resolve_for_this_platform(monkeypatch, tmp_path):
    """Read target-scoped preparation the way Pixi resolves it for the host."""
    manifest = tmp_path / "pyproject.toml"
    manifest.write_text(
        "[tool.pixi.environments]\n"
        'scoped = { features = ["scoped"], no-default-feature = true }\n'
        "inherited = { features = [] }\n"
        "[tool.pixi.tasks]\n"
        'prepare = "echo workspace"\n'
        "[tool.pixi.target.linux-64.tasks]\n"
        'prepare = "echo workspace-linux"\n'
        "[tool.pixi.feature.scoped.tasks]\n"
        'prepare = "echo feature"\n'
        "[tool.pixi.feature.scoped.target.linux-64.tasks]\n"
        'prepare = "echo feature-linux"\n'
    )
    monkeypatch.setattr(
        environments, "settings", Settings(environment_manifest=manifest)
    )
    monkeypatch.setattr(environments.platform, "system", lambda: "Linux")
    monkeypatch.setattr(environments.platform, "machine", lambda: "x86_64")
    assert environments.prepare_environments() == ("scoped", "inherited")
    # A target table overrides the unscoped declaration in the same table.
    assert environments.prepare_command("scoped") == json.dumps(["echo feature-linux"])
    assert environments.prepare_command("inherited") == json.dumps(
        ["echo workspace-linux"]
    )
    # Another platform's target table contributes nothing here.
    monkeypatch.setattr(environments.platform, "machine", lambda: "arm64")
    assert environments.prepare_command("scoped") == json.dumps(["echo feature"])
    assert environments.prepare_command("inherited") == json.dumps(["echo workspace"])
    # An unknown platform keeps the unscoped declarations instead of guessing.
    monkeypatch.setattr(environments.platform, "machine", lambda: "riscv64")
    assert environments.prepare_command("inherited") == json.dumps(["echo workspace"])


@pytest.mark.parametrize("name", ["../escape", "/absolute", "-option", ""])
def test_environment_names_cannot_escape_preparation_directory(
    monkeypatch, tmp_path, name
):
    """Reject names that could write outside the manifest's preparation directory."""
    manifest = tmp_path / "pyproject.toml"
    manifest.write_text(f"[tool.pixi.environments]\n{json.dumps(name)} = []\n")
    monkeypatch.setattr(
        environments, "settings", Settings(environment_manifest=manifest)
    )
    with pytest.raises(ValueError, match="Invalid Pixi environment name"):
        environments.declared_environments()
