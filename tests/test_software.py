"""Package-source delegation, locked environments, and native argument preservation."""

import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from click import unstyle
from typer.testing import CliRunner

from biov import environments, software
from biov.cli import app
from biov.config import Settings


@pytest.fixture(autouse=True)
def managed_pixi_only(monkeypatch, tmp_path):
    """Resolve the managed copy unless a test provides its own Pixi candidate."""
    monkeypatch.delenv(environments.PIXI_BIN, raising=False)
    monkeypatch.setenv("BIOV_ENVIRONMENT_ROOT", str(tmp_path / "manager"))
    monkeypatch.setattr(environments.shutil, "which", lambda name: None)
    monkeypatch.setattr(
        environments, "settings", Settings(environment_root=tmp_path / "manager")
    )


def write_manifest(directory: Path, *, prepare: str | None = None) -> Path:
    """Write a project manifest declaring two environments, plus its lock.

    Returns:
        Manifest path a Pixi selection points at.
    """
    tasks = (
        f'[tool.pixi.feature.mafft.tasks]\nprepare = "{prepare}"\n'
        if prepare is not None
        else ""
    )
    manifest = directory / "pyproject.toml"
    manifest.write_text(
        "[tool.pixi.environments]\n"
        'mafft = { features = ["mafft"], no-default-feature = true }\n'
        "plain = { features = [], no-default-feature = true }\n"
        f"{tasks}"
    )
    manifest.with_name("pixi.lock").write_text("lock = 1\n")
    return manifest


def preparation_records(root: Path) -> list[Path]:
    """Return every preparation record below one directory.

    Returns:
        Record paths, wherever Pixi put the environments they describe.
    """
    return list(root.rglob(environments.PREPARED_RECORD))


def fake_pixi_run(argv: list[str], returncode: int = 0) -> subprocess.CompletedProcess:
    """Answer one faked Pixi invocation as the pinned manager would.

    Returns:
        Synthetic process result; an install creates its environment prefix and
        an info call reports the prefixes that exist under the project.
    """
    if "--manifest-path" not in argv:
        return subprocess.CompletedProcess(argv, returncode)
    manifest = Path(argv[argv.index("--manifest-path") + 1])
    root = manifest.parent / ".pixi" / "envs"
    if argv[1] == "install":
        (root / argv[argv.index("--environment") + 1]).mkdir(
            parents=True, exist_ok=True
        )
    if argv[1] == "info":
        listing = [
            {"name": prefix.name, "prefix": str(prefix)}
            for prefix in sorted(root.glob("*"))
            if prefix.is_dir()
        ]
        return subprocess.CompletedProcess(
            argv, returncode, json.dumps({"environments_info": listing})
        )
    return subprocess.CompletedProcess(argv, returncode)


def record_pixi_runs(monkeypatch, returncode: int = 0) -> list[list[str]]:
    """Record every Pixi argv and report success.

    Returns:
        The list the recorded argv lists are appended to.
    """
    calls: list[list[str]] = []
    executable = environments.managed_pixi_path()
    executable.parent.mkdir(parents=True, exist_ok=True)
    executable.touch()

    def run(argv, **kwargs):
        calls.append(list(argv))
        return fake_pixi_run(argv, returncode)

    monkeypatch.setattr(subprocess, "run", run)
    return calls


@pytest.mark.parametrize("source", ["override", "path"])
def test_relative_pixi_path_survives_execution_cwd(monkeypatch, tmp_path, source):
    """Resolve the manager in the caller's directory before changing cwd."""
    monkeypatch.chdir(tmp_path)
    executable = tmp_path / "manager" / "pixi"
    executable.parent.mkdir()
    executable.touch()
    working_directory = tmp_path / "analysis"
    working_directory.mkdir()
    manifest = write_manifest(tmp_path)
    config = Settings(environment_manifest=manifest)
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    if source == "override":
        monkeypatch.setenv(environments.PIXI_BIN, "manager/pixi")
    else:
        monkeypatch.setattr(environments.shutil, "which", lambda name: "manager/pixi")

    def run(argv, **kwargs):
        assert argv[0] == str(executable)
        if argv[1:] == ["--version"]:
            return subprocess.CompletedProcess(
                argv, 0, f"pixi {environments.PIXI_VERSION}"
            )
        if argv[1] == "install":
            return subprocess.CompletedProcess(argv, 0)
        assert kwargs["cwd"] == working_directory
        return subprocess.CompletedProcess(argv, 23)

    monkeypatch.setattr(subprocess, "run", run)
    result = CliRunner().invoke(
        app, ["exec", "--cwd", str(working_directory), "mafft", "--version"]
    )
    assert result.exit_code == 23, result.output


@pytest.mark.parametrize("tool", ["mafft", "conda:mafft"])
def test_pixi_routing_arguments_and_inherited_environment(monkeypatch, tmp_path, tool):
    """Keep native arguments and caller caches while selecting the declared env."""
    manifest = write_manifest(tmp_path)
    config = Settings.model_validate(
        {
            "home": tmp_path / "cache",
            "environment_root": tmp_path / "environments",
            "environment_manifest": manifest,
        }
    )
    monkeypatch.setattr(software, "settings", config)
    monkeypatch.setattr("biov.environments.settings", config)
    executable = environments.managed_pixi_path()
    executable.parent.mkdir(parents=True)
    executable.touch()
    for key in ("FSSPEC_FILECACHE", "HF_HOME", "TORCH_HOME", "XDG_CACHE_HOME"):
        monkeypatch.setenv(key, "caller-setting")
    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 23)

    monkeypatch.setattr(subprocess, "run", run)
    arguments = ("analysis.py", "a;$(touch forbidden)", "one\ntwo", "")
    result = CliRunner().invoke(app, ["exec", "--cwd", str(tmp_path), tool, *arguments])
    assert result.exit_code == 23, result.output
    argv, kwargs = calls.pop()
    assert argv == [
        str(environments.pixi_path()),
        "run",
        "--no-config",
        "--frozen",
        "--manifest-path",
        str(manifest),
        "--environment",
        "mafft",
        "--",
        "mafft",
        *arguments,
    ]
    assert calls.pop()[0] == [
        str(environments.pixi_path()),
        "install",
        "--no-config",
        "--manifest-path",
        str(manifest),
        "--locked",
        "--environment",
        "mafft",
    ]
    assert kwargs["cwd"] == tmp_path.resolve()
    assert kwargs["env"] == {
        **os.environ,
        "BIOV_HOME": str(tmp_path / "cache"),
        "BIOV_CACHE_HTTP": "true",
    }
    assert kwargs["check"] is False
    assert not config.home.exists()
    software.run_software("unlisted", ("--help",), cwd=tmp_path)
    assert calls.pop()[0] == [
        str(executable),
        "exec",
        "-s",
        "unlisted",
        "--",
        "unlisted",
        "--help",
    ]
    # Only a name declared by the selected manifest owns an environment.
    software.run_software("samtools", ("--version",), cwd=tmp_path)
    assert calls.pop()[0] == [
        str(executable),
        "exec",
        "-s",
        "samtools",
        "--",
        "samtools",
        "--version",
    ]
    assert not calls


@pytest.mark.parametrize("entry", ["workspace-task", "feature-task", "executable"])
def test_declared_entry_preserves_quotes_through_native_task_parsing(
    monkeypatch, tmp_path, entry
):
    """Apply Pixi task quoting without changing executable argv."""
    manifest = write_manifest(tmp_path)
    config = Settings(environment_manifest=manifest)
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    calls = record_pixi_runs(monkeypatch)
    task = {"name": "mafft", "cmd": "entry-point"}
    listing = [
        {"environment": "plain", "tasks": [task], "features": []},
        {
            "environment": "mafft",
            "tasks": [task] if entry == "workspace-task" else [],
            "features": [
                {
                    "name": "mafft",
                    "tasks": [task] if entry == "feature-task" else [],
                }
            ],
        },
    ]
    arguments = [
        "single'quote",
        'double"quote',
        "paired'quotes'",
        '"paired quotes"',
        "literal;'$(touch forbidden)'",
        "one\ntwo",
        "",
        "--help",
        "trailing\\",
    ]

    def run(argv, **kwargs):
        calls.append(list(argv))
        if argv[1] == "task":
            return subprocess.CompletedProcess(argv, 0, json.dumps(listing).encode())
        assert argv[1] == "run"
        assert kwargs["cwd"] == tmp_path
        received = argv[argv.index("--") + 1 :]
        if entry != "executable":
            assert "--executable" not in argv
            # Pixi 0.81.0 as_script surrounds each extra argument with single
            # quotes; parse that script rather than mirroring BioV's escaping.
            received = shlex.split(
                "mafft " + " ".join(f"'{arg}'" for arg in received[1:])
            )
        assert received == ["mafft", *arguments]
        return subprocess.CompletedProcess(argv, 37)

    monkeypatch.setattr(subprocess, "run", run)
    result = CliRunner().invoke(
        app,
        ["exec", "--no-install", "--cwd", str(tmp_path), "mafft", *arguments],
    )
    assert result.exit_code == 37, result.output
    assert [argv[1] for argv in calls] == ["task", "run"]


@pytest.mark.parametrize(
    ("name", "entry"),
    [
        ("mafft", "missing"),
        ("plain", "missing"),
        ("mafft", "task"),
        ("mafft", "executable"),
    ],
)
def test_declared_missing_entry_is_actionable_without_rewriting_real_exit_127(
    monkeypatch, tmp_path, name, entry
):
    """Explain absent entries and preserve a real task or executable's status."""
    manifest = write_manifest(tmp_path)
    config = Settings(environment_manifest=manifest)
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    calls = record_pixi_runs(monkeypatch)
    activated_path = str(tmp_path / "activated" / "bin")
    looked_up = []

    def lookup(command, *, path=None):
        if path is None:
            return None
        assert command == name and path == activated_path
        looked_up.append(command)
        return str(Path(path) / name) if entry == "executable" else None

    def run(argv, **kwargs):
        calls.append(list(argv))
        if argv[1] == "run":
            return subprocess.CompletedProcess(argv, 127)
        if argv[1] == "task":
            listing = [
                {
                    "environment": name,
                    "tasks": [],
                    "features": [
                        {
                            "name": name,
                            "tasks": [{"name": name}] if entry == "task" else [],
                        }
                    ],
                }
            ]
            return subprocess.CompletedProcess(argv, 0, json.dumps(listing).encode())
        assert argv[1] == "shell-hook"
        assert "--as-is" in argv
        activation = {
            "environment_variables": {"PATH": activated_path},
            "activation_scripts": [],
        }
        return subprocess.CompletedProcess(argv, 0, json.dumps(activation).encode())

    monkeypatch.setattr(software.shutil, "which", lookup)
    monkeypatch.setattr(subprocess, "run", run)
    result = CliRunner().invoke(app, ["exec", "--no-install", name])
    if entry == "missing":
        assert result.exit_code == 2, result.output
        assert (
            f"declared environment '{name}' declares no '{name}' task or executable"
            in result.output
        )
        assert "add one" in result.output and "biov setup" in result.output
        assert (
            'feature."mafft".tasks' if name == "mafft" else "select a feature"
        ) in result.output
    else:
        assert result.exit_code == 127, result.output
        assert "declares no" not in result.output
    assert [argv[1] for argv in calls] == (
        ["run", "task"] if entry == "task" else ["run", "task", "shell-hook"]
    )
    assert looked_up == ([] if entry == "task" else [name])
    assert not preparation_records(tmp_path)


@pytest.mark.parametrize("unavailable", ["task-list", "activation"])
def test_unavailable_pixi_introspection_preserves_native_exit_status(
    monkeypatch, tmp_path, unavailable
):
    """Keep an entry's own status when Pixi cannot report tasks or activation."""
    manifest = write_manifest(tmp_path)
    config = Settings(environment_manifest=manifest)
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    executable = environments.managed_pixi_path()
    executable.parent.mkdir(parents=True, exist_ok=True)
    executable.touch()
    calls: list[str] = []

    def run(argv, **kwargs):
        calls.append(argv[1])
        if argv[1] == "run":
            return subprocess.CompletedProcess(argv, 127)
        if argv[1] == "task" and unavailable == "activation":
            return subprocess.CompletedProcess(argv, 0, json.dumps([]).encode())
        raise subprocess.CalledProcessError(1, argv)

    monkeypatch.setattr(subprocess, "run", run)
    result = CliRunner().invoke(app, ["exec", "--no-install", "mafft"])

    assert result.exit_code == 127, result.output
    assert "declares no" not in result.output
    assert calls == (
        ["run", "task"] if unavailable == "task-list" else ["run", "task", "shell-hook"]
    )


def test_missing_runtime_and_directory_fail_without_fallback(monkeypatch, tmp_path):
    """A missing runtime fails once; invalid cwd fails before process creation."""
    config = Settings.model_validate(
        {
            "environment_manifest": tmp_path / "pyproject.toml",
        }
    )
    write_manifest(tmp_path)
    monkeypatch.setattr(software, "settings", config)
    monkeypatch.setattr("biov.environments.settings", config)
    calls = []

    def run(argv, **kwargs):
        calls.append(argv)
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(subprocess, "run", run)
    with pytest.raises(FileNotFoundError, match="biov setup"):
        software.run_software("mafft", ("--version",), cwd=tmp_path)
    assert not calls
    with pytest.raises(FileNotFoundError):
        software.run_software("mafft", cwd=tmp_path / "missing")
    file = tmp_path / "file"
    file.touch()
    with pytest.raises(NotADirectoryError):
        software.run_software("mafft", cwd=file)
    assert not calls


@pytest.mark.parametrize(
    ("tool", "expected"),
    [
        ("seqtk", ["exec", "-s", "seqtk", "--", "seqtk"]),
        ("conda:seqtk", ["exec", "-s", "seqtk", "--", "seqtk"]),
        ("pypi:ruff", ["tool", "run", "ruff"]),
        ("npm:@scope/tool", ["--yes", "@scope/tool"]),
    ],
)
def test_package_sources_preserve_arguments_caches_cwd_and_status(
    monkeypatch, tmp_path, tool, expected
):
    """Delegate package execution without resolving another source or a workspace."""
    manifest = write_manifest(tmp_path)
    config = Settings(
        home=tmp_path / "cache",
        cache_http=False,
        environment_root=tmp_path / "environments",
        environment_manifest=manifest,
    )
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    monkeypatch.setattr(
        software.shutil,
        "which",
        lambda name: str(tmp_path / name) if name in {"uv", "npx"} else None,
    )
    executable = environments.managed_pixi_path()
    executable.parent.mkdir(parents=True)
    executable.touch()
    for key in ("FSSPEC_FILECACHE", "HF_HOME", "TORCH_HOME", "XDG_CACHE_HOME"):
        monkeypatch.setenv(key, "caller-setting")
    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 23)

    monkeypatch.setattr(subprocess, "run", run)
    arguments = ["--help", "a;$(touch forbidden)", "one\ntwo", ""]
    result = CliRunner().invoke(app, ["exec", "--cwd", str(tmp_path), tool, *arguments])
    assert result.exit_code == 23, result.output
    manager = (
        str(environments.managed_pixi_path())
        if not tool.startswith(("pypi:", "npm:"))
        else str(tmp_path / ("uv" if tool.startswith("pypi:") else "npx"))
    )
    assert calls == [
        (
            [manager, *expected, *arguments],
            {
                "cwd": tmp_path,
                "env": {
                    **os.environ,
                    "BIOV_HOME": str(config.home),
                    "BIOV_CACHE_HTTP": "false",
                },
                "check": False,
            },
        )
    ]
    assert not config.home.exists()
    assert not preparation_records(tmp_path)


@pytest.mark.parametrize(("source", "manager"), [("pypi", "uv"), ("npm", "npx")])
def test_missing_source_manager_reports_installation_without_fallback(
    monkeypatch, tmp_path, source, manager
):
    """An absent source manager reports the remedy without starting a process."""
    config = Settings(environment_manifest=write_manifest(tmp_path))
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    calls = record_pixi_runs(monkeypatch)
    with pytest.raises(FileNotFoundError, match=manager) as failure:
        software.run_software(f"{source}:mafft", ("--version",), cwd=tmp_path)
    assert "install" in str(failure.value).lower()
    result = CliRunner().invoke(app, ["exec", f"{source}:mafft", "--version"])
    assert result.exit_code == 2, result.output
    assert manager in result.output and "install" in result.output.lower()
    assert "Traceback" not in result.output
    assert not calls


@pytest.mark.parametrize(
    "tool", ["unlisted", "conda:unlisted", "pypi:unlisted", "npm:unlisted"]
)
def test_transient_sources_reject_no_install_before_manager_lookup(
    monkeypatch, tmp_path, tool
):
    """Never invoke a package runner when installation has been forbidden."""
    config = Settings(environment_manifest=write_manifest(tmp_path))
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)

    def lookup(name):
        pytest.fail(f"Must reject --no-install before looking up {name}")

    monkeypatch.setattr(software.shutil, "which", lookup)
    calls = record_pixi_runs(monkeypatch)
    with pytest.raises(ValueError, match="no-install"):
        software.run_software(tool, cwd=tmp_path, install=False)
    result = CliRunner().invoke(app, ["exec", "--no-install", tool])
    assert result.exit_code == 2, result.output
    assert "no-install" in result.output and not calls


@pytest.mark.parametrize("source", ["conda", "pypi", "npm"])
@pytest.mark.parametrize("name", ["", "--help"])
def test_recognized_source_requires_a_package_coordinate(
    monkeypatch, tmp_path, source, name
):
    """Reject missing coordinates and source-manager options at the boundary."""
    monkeypatch.setattr(software, "settings", Settings(home=tmp_path / "cache"))
    calls = record_pixi_runs(monkeypatch)
    with pytest.raises(ValueError):
        software.run_software(f"{source}:{name}", cwd=tmp_path)
    result = CliRunner().invoke(app, ["exec", f"{source}:{name}"])
    assert result.exit_code == 2, result.output
    assert not calls


@pytest.mark.parametrize(
    "tool",
    [
        "other:tool",
        "CONDA:tool",
        "https://example.test/tool",
        "conda",
        "pypi",
        "npm",
        "echo",
    ],
)
def test_unrecognized_coordinate_uses_conda_without_host_fallback(
    monkeypatch, tmp_path, tool
):
    """Keep the whole coordinate and use Conda even when a host command exists."""
    config = Settings(environment_manifest=write_manifest(tmp_path))
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    monkeypatch.setattr(
        software.shutil, "which", lambda name: "/bin/echo" if name == "echo" else None
    )
    calls = record_pixi_runs(monkeypatch, returncode=19)
    result = CliRunner().invoke(app, ["exec", tool, "--version"])
    assert result.exit_code == 19, result.output
    assert calls == [
        [
            str(environments.managed_pixi_path()),
            "exec",
            "-s",
            tool,
            "--",
            tool,
            "--version",
        ]
    ]


@pytest.mark.parametrize("delimiter", [[], ["--"]])
def test_exec_passes_native_options_after_coordinate(monkeypatch, tmp_path, delimiter):
    """Program options keep their meaning even when BioV has the same option."""
    monkeypatch.setattr(software, "settings", Settings(home=tmp_path / "cache"))
    calls = record_pixi_runs(monkeypatch, returncode=17)
    arguments = ["--help", "--cwd", "native-work", "--no-install", "--", ""]
    result = CliRunner().invoke(
        app, ["exec", "--cwd", str(tmp_path), "conda:echo", *delimiter, *arguments]
    )
    assert result.exit_code == 17, result.output
    assert calls == [
        [
            str(environments.managed_pixi_path()),
            "exec",
            "-s",
            "echo",
            "--",
            "echo",
            *arguments,
        ]
    ]


def test_pixi_passthrough_forwards_arguments_and_exit_status(monkeypatch, tmp_path):
    """Run the resolved manager with native arguments, and never install one."""
    config = Settings(environment_root=tmp_path / "environments")
    monkeypatch.setattr(environments, "settings", config)
    result = CliRunner().invoke(
        app,
        ["pixi", "run", "--environment", "samtools", "--", "samtools", "--version"],
    )
    assert result.exit_code == 2, result.output
    assert "biov setup" in result.output
    executable = environments.managed_pixi_path()
    executable.parent.mkdir(parents=True)
    executable.touch()
    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 7)

    monkeypatch.setattr(subprocess, "run", run)
    result = CliRunner().invoke(
        app,
        ["pixi", "run", "--environment", "samtools", "--", "samtools", "--version"],
    )
    assert result.exit_code == 7, result.output
    assert calls == [
        (
            [
                str(executable),
                "run",
                "--environment",
                "samtools",
                "--",
                "samtools",
                "--version",
            ],
            {"check": False},
        )
    ]


def test_pixi_cli_passes_native_help_and_delimiter_to_process(monkeypatch, tmp_path):
    """Pass exact argv, including bare help, to a real Pixi-like process."""
    recorded = tmp_path / "pixi-argv.json"
    script = tmp_path / "fake-pixi.py"
    script.write_text(
        "import json\n"
        "import sys\n"
        "from pathlib import Path\n"
        f"Path({str(recorded)!r}).write_text(json.dumps(sys.argv[1:]))\n"
        "raise SystemExit(17)\n"
    )
    monkeypatch.setattr("biov.cli.pixi_command", lambda: (sys.executable, str(script)))
    for native_arguments in (
        ["run", "samtools", "--help"],
        ["run", "--", "echo", "--help"],
        ["run", "--environment", "samtools", "--", "samtools", "--version"],
        ["--", "--help"],
        ["--help"],
    ):
        recorded.write_text("not run")
        result = CliRunner().invoke(app, ["pixi", *native_arguments])
        assert result.exit_code == 17, result.output
        assert json.loads(recorded.read_text()) == native_arguments
    recorded.write_text("not run")
    help_result = CliRunner().invoke(app, ["--help"])
    assert help_result.exit_code == 0, help_result.output
    assert "pixi" in help_result.output
    assert recorded.read_text() == "not run"


def test_execution_prefers_an_installed_pixi(monkeypatch, tmp_path):
    """Run a declared environment through a PATH Pixi of the pinned version."""
    manifest = write_manifest(tmp_path)
    config = Settings(environment_manifest=manifest)
    monkeypatch.setattr(software, "settings", config)
    monkeypatch.setattr(environments, "settings", config)
    system = tmp_path / "bin" / "pixi"
    system.parent.mkdir()
    system.write_bytes(b"executable")
    monkeypatch.setattr(environments.shutil, "which", lambda name: str(system))
    calls = []

    def run(argv, **kwargs):
        calls.append(argv)
        return subprocess.CompletedProcess(
            argv, 0, f"pixi {environments.PIXI_VERSION}\n"
        )

    monkeypatch.setattr(subprocess, "run", run)
    result = software.run_software("mafft", ("--version",), cwd=tmp_path)
    assert result.returncode == 0
    assert calls == [
        [str(system), "--version"],
        [
            str(system),
            "install",
            "--no-config",
            "--manifest-path",
            str(manifest),
            "--locked",
            "--environment",
            "mafft",
        ],
        [
            str(system),
            "run",
            "--no-config",
            "--frozen",
            "--manifest-path",
            str(manifest),
            "--environment",
            "mafft",
            "--",
            "mafft",
            "--version",
        ],
    ]


def test_shipped_environment_provisions_instead_of_command_not_found(
    monkeypatch, tmp_path
):
    """Install the declared environment instead of reporting a missing command."""
    config = Settings(environment_root=tmp_path / "environments")
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    assert "mafft" in environments.declared_environments()
    calls = record_pixi_runs(monkeypatch)
    result = CliRunner().invoke(app, ["exec", "mafft", "--", "--version"])
    assert result.exit_code == 0, result.output
    executable = str(environments.pixi_path())
    assert calls == [
        [
            executable,
            "install",
            "--no-config",
            "--manifest-path",
            str(environments.environment_manifest()),
            "--locked",
            "--environment",
            "mafft",
        ],
        [
            executable,
            "run",
            "--no-config",
            "--frozen",
            "--manifest-path",
            str(environments.environment_manifest()),
            "--environment",
            "mafft",
            "--",
            "mafft",
            "--version",
        ],
    ]
    # The bundled manifest and lock are published before Pixi reads them, so a
    # fresh device needs no separate setup step.
    manifest = environments.environment_manifest()
    assert manifest.is_file() and manifest.with_name("pixi.lock").is_file()
    # mafft is never executed as a bare host command, which is what used to hide
    # "this environment was never installed" behind "command not found".
    assert not any(argv == ["mafft", "--version"] for argv in calls)


def test_exec_provisions_and_prepares_once(monkeypatch, tmp_path):
    """Install from the lock, prepare once, then run the command itself."""
    manifest = write_manifest(tmp_path, prepare="git fetch uce")
    config = Settings.model_validate(
        {
            "environment_root": tmp_path / "environments",
            "environment_manifest": manifest,
        }
    )
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    calls = record_pixi_runs(monkeypatch)
    executable = str(environments.pixi_path())
    install = ["install", "--no-config", "--manifest-path", str(manifest), "--locked"]
    run = [
        "run",
        "--no-config",
        "--frozen",
        "--manifest-path",
        str(manifest),
    ]
    software.run_software("mafft", ("--version",), cwd=tmp_path)
    assert calls == [
        [executable, *install, "--environment", "mafft"],
        [
            executable,
            "info",
            "--no-config",
            "--json",
            "--manifest-path",
            str(manifest),
        ],
        [
            executable,
            "run",
            "--no-config",
            "--manifest-path",
            str(manifest),
            "--frozen",
            "--environment",
            "mafft",
            "prepare",
        ],
        [executable, *run, "--environment", "mafft", "--", "mafft", "--version"],
    ]
    records = preparation_records(tmp_path)
    assert len(records) == 1
    assert "environment = mafft" in records[0].read_text()
    # Preparation takes minutes, so a second exec installs and runs but does not
    # repeat it.
    calls.clear()
    software.run_software("mafft", ("--version",), cwd=tmp_path)
    assert [argv[1] for argv in calls] == ["install", "info", "run"]
    # Editing the prepare task invalidates the record.
    write_manifest(tmp_path, prepare="git fetch another")
    calls.clear()
    software.run_software("mafft", ("--version",), cwd=tmp_path)
    assert [argv[1] for argv in calls] == ["install", "info", "run", "run"]
    assert calls[2][-1] == "prepare"
    # So does a new lock.
    manifest.with_name("pixi.lock").write_text("lock = 2\n")
    calls.clear()
    software.run_software("mafft", ("--version",), cwd=tmp_path)
    assert [argv[1] for argv in calls] == ["install", "info", "run", "run"]
    # Pixi cleaning the environment takes the record with the prefix, so the
    # reinstalled environment is prepared instead of reused as ready.
    shutil.rmtree(manifest.parent / ".pixi" / "envs" / "mafft")
    calls.clear()
    software.run_software("mafft", ("--version",), cwd=tmp_path)
    assert [argv[1] for argv in calls] == ["install", "info", "run", "run"]
    assert calls[2][-1] == "prepare"


def test_failed_prepare_records_no_marker(monkeypatch, tmp_path):
    """A preparation that fails is retried instead of hidden by a marker."""
    manifest = write_manifest(tmp_path, prepare="git fetch uce")
    config = Settings.model_validate(
        {
            "environment_root": tmp_path / "environments",
            "environment_manifest": manifest,
        }
    )
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    executable = environments.managed_pixi_path()
    executable.parent.mkdir(parents=True)
    executable.touch()
    calls = []

    def run(argv, **kwargs):
        calls.append(list(argv))
        if argv[1] == "run" and argv[-1] == "prepare":
            raise subprocess.CalledProcessError(1, argv)
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(subprocess, "run", run)
    with pytest.raises(subprocess.CalledProcessError):
        software.run_software("mafft", ("--version",), cwd=tmp_path)
    assert [argv[1] for argv in calls] == ["install", "info", "run"]
    assert not preparation_records(tmp_path)
    # The failed preparation is retried instead of being treated as done.
    calls.clear()

    def succeeded(argv, **kwargs):
        calls.append(list(argv))
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(subprocess, "run", succeeded)
    software.run_software("mafft", ("--version",), cwd=tmp_path)
    assert [argv[1] for argv in calls] == ["install", "info", "run", "run"]


@pytest.mark.parametrize("tool", ["mafft", "conda:mafft"])
@pytest.mark.parametrize("force_color", [False, True])
def test_exec_no_install_never_touches_the_environment(
    monkeypatch, tmp_path, tool, force_color
):
    """--no-install skips installation and preparation and runs with --as-is."""
    manifest = write_manifest(tmp_path, prepare="git fetch uce")
    config = Settings.model_validate(
        {
            "environment_root": tmp_path / "environments",
            "environment_manifest": manifest,
        }
    )
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    calls = record_pixi_runs(monkeypatch, returncode=29)
    arguments = ["exec", "--no-install", tool, "--version"]
    result = CliRunner().invoke(app, arguments)
    assert result.exit_code == 29, result.output
    assert calls == [
        [
            str(environments.pixi_path()),
            "run",
            "--no-config",
            "--as-is",
            "--manifest-path",
            str(manifest),
            "--environment",
            "mafft",
            "--",
            "mafft",
            "--version",
        ]
    ]
    assert not preparation_records(tmp_path)
    help_result = CliRunner().invoke(
        app,
        ["exec", "--help"],
        env={
            "FORCE_COLOR": "1" if force_color else "",
            "NO_COLOR": "" if force_color else "1",
            "TERM": "xterm" if force_color else "dumb",
        },
    )
    assert help_result.exit_code == 0, help_result.output
    if force_color:
        assert "\x1b[" in help_result.output
    assert "--no-install" in unstyle(help_result.output)


def test_no_install_missing_bundled_workspace_reports_setup(monkeypatch, tmp_path):
    """Explain a missing bundled workspace without creating or installing it."""
    config = Settings(environment_root=tmp_path / "environments")
    monkeypatch.setattr(environments, "settings", config)
    monkeypatch.setattr(software, "settings", config)
    calls = record_pixi_runs(monkeypatch)
    result = CliRunner().invoke(app, ["exec", "--no-install", "r", "--help"])
    assert result.exit_code == 2, result.output
    assert "no initialized workspace" in result.output
    assert "biov setup r" in result.output
    assert not (config.environment_root / "workspaces").exists()
    assert not calls


def test_no_install_is_forwarded_to_the_execution_host(monkeypatch, tmp_path):
    """The escape hatch reaches the host that owns the environment."""
    config = Settings(execution_host="compute", home=tmp_path / "cache")
    monkeypatch.setattr(software, "settings", config)
    calls = record_pixi_runs(monkeypatch)
    software.run_software("mafft", ("--version",), install=False)
    assert shlex.split(calls[0][-1]) == [
        "biov",
        "exec",
        "--no-install",
        "--",
        "mafft",
        "--version",
    ]
    calls.clear()
    software.run_software("mafft", ("--version",))
    assert shlex.split(calls[0][-1]) == ["biov", "exec", "--", "mafft", "--version"]
