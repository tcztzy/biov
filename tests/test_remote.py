"""Package execution settings and transparent SSH routing without remote connections."""

import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from typer.testing import CliRunner

from biov import config as config_module
from biov import software
from biov.cli import app
from biov.config import Settings
from biov.remote import run_remote


@pytest.fixture
def fake_uv(monkeypatch, tmp_path):
    """Execute a Python program through a real process without fetching a package."""
    executable = tmp_path / "bin" / "uv"
    executable.parent.mkdir()
    executable.write_text(
        f"#!{sys.executable}\n"
        "import os,sys\n"
        'assert sys.argv[1:4] == ["tool", "run", "python"]\n'
        "os.execv(sys.executable, [sys.executable, *sys.argv[4:]])\n"
    )
    executable.chmod(0o700)
    monkeypatch.setenv("PATH", f"{executable.parent}{os.pathsep}{os.environ['PATH']}")


def test_exec_reaches_execution_host_with_private_application_settings(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fake_uv
) -> None:
    """Simulate SSH transport while actually running the receiving BioV process."""
    remote_config = tmp_path / "remote.toml"
    remote_config.write_text(f"home = {json.dumps(str(tmp_path / 'remote-cache'))}\n")
    selected = Settings.model_validate(
        {
            "execution_host": "compute",
            "execution_cwd": tmp_path,
            "home": tmp_path / "local-only-cache",
        }
    )
    monkeypatch.setattr(software, "settings", selected)
    native_run = subprocess.run

    def ssh(argv, **kwargs):
        assert argv[:4] == ["ssh", "-T", "--", "compute"]
        return native_run(
            ["/bin/sh", "-c", argv[-1]],
            cwd=tmp_path,
            env={
                **{k: v for k, v in os.environ.items() if not k.startswith("BIOV_")},
                "PATH": f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}",
                "BIOV_CONFIG": str(remote_config),
            },
            check=False,
        )

    monkeypatch.setattr(subprocess, "run", ssh)
    result = CliRunner().invoke(
        app,
        [
            "exec",
            "--",
            "pypi:python",
            "-c",
            "from pathlib import Path; import os,sys; Path('remote-result.txt').write_text(os.environ['BIOV_HOME'] + '\\n' + sys.argv[1]); sys.exit(29)",
            "$(touch forbidden)",
        ],
    )
    assert result.exit_code == 29, result.output
    assert (
        tmp_path / "remote-result.txt"
    ).read_text() == f"{tmp_path}/remote-cache\n$(touch forbidden)"
    assert not (tmp_path / "forbidden").exists()


def test_application_config_and_environment_precedence(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Read TOML application settings and permit operator environment overrides."""
    config = tmp_path / "config.toml"
    config.write_text(
        'execution_host = "compute"\nexecution_cwd = "/remote/work"\n'
        'environment_root = "/remote/environments"\n'
    )
    monkeypatch.setenv("BIOV_CONFIG", str(config))
    monkeypatch.setenv("BIOV_CACHE_HTTP", "false")
    settings = Settings()
    assert settings.execution_host == "compute"
    assert settings.execution_cwd == Path("/remote/work")
    assert settings.environment_root == Path("/remote/environments")
    assert settings.cache_http is False
    monkeypatch.setenv("BIOV_EXECUTION_HOST", "another-host")
    assert Settings().execution_host == "another-host"
    monkeypatch.delenv("BIOV_CONFIG")
    monkeypatch.setattr("biov.config.user_config_path", lambda _: tmp_path)
    assert Settings().environment_root == Path("/remote/environments")


@pytest.mark.parametrize("option", ["--config", "--config="])
@pytest.mark.parametrize("selection", ["missing", "directory", "empty"])
def test_cli_config_overrides_environment_before_import(
    tmp_path, option, selection, fake_uv
):
    """Load the root CLI selection before an unusable environment selection."""
    config = tmp_path / "config.toml"
    cache = tmp_path / "cache"
    config.write_text(f"home = {json.dumps(str(cache))}\n")
    invalid = {
        "missing": str(tmp_path / "missing"),
        "directory": str(tmp_path),
        "empty": "",
    }
    biov = shutil.which("biov", path=str(Path(sys.executable).parent))
    assert biov is not None
    arguments = [option, str(config)] if option == "--config" else [f"{option}{config}"]
    native_config = str(tmp_path / "native-config")
    result = subprocess.run(
        [
            biov,
            *arguments,
            "exec",
            "--",
            "pypi:python",
            "-c",
            "import json,os,sys; print(json.dumps([os.environ['BIOV_HOME'], sys.argv[1:]]))",
            "--config",
            native_config,
        ],
        env={
            **{k: v for k, v in os.environ.items() if not k.startswith("BIOV_")},
            "BIOV_CONFIG": invalid[selection],
        },
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == [str(cache), ["--config", native_config]]


def test_unusable_config_still_fails_for_cli_and_python_import(tmp_path):
    """Reject invalid selections without interpreting a Python program's argv."""
    config = tmp_path / "config.toml"
    config.touch()
    biov = shutil.which("biov", path=str(Path(sys.executable).parent))
    assert biov is not None
    for command, status in (
        ([biov, "--help"], 2),
        ([sys.executable, "-c", "import biov", "--config", str(config)], 1),
    ):
        result = subprocess.run(
            command,
            env={**os.environ, "BIOV_CONFIG": str(tmp_path / "missing")},
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == status
        assert "BIOV_CONFIG" in result.stderr and "--config" in result.stderr
        if status == 2:
            assert "Traceback" not in result.stderr


def test_select_config_file_exports_the_selection(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """An explicit selection becomes the file child processes inherit."""
    config = tmp_path / "config.toml"
    config.write_text("cache_http = false\n")
    monkeypatch.setenv("BIOV_CONFIG", str(tmp_path / "stale.toml"))
    monkeypatch.setattr(
        config_module.settings, "__dict__", dict(vars(config_module.settings))
    )

    config_module.select_config_file(config)

    assert os.environ["BIOV_CONFIG"] == str(config.resolve())
    assert config_module.settings.cache_http is False


def test_cli_config_reaches_an_executed_analysis_script(tmp_path: Path) -> None:
    """Run a script that imports BioV under its own invocation's selection."""
    config = tmp_path / "selected.toml"
    cache = tmp_path / "selected-cache"
    config.write_text(f"home = {json.dumps(str(cache))}\n")
    script = tmp_path / "settings_report.py"
    script.write_text(
        "import json\n"
        "from biov.config import settings\n"
        "print(json.dumps({'config': str(settings.config),"
        " 'home': str(settings.home)}))\n"
    )
    biov = shutil.which("biov", path=str(Path(sys.executable).parent))
    assert biov is not None

    result = subprocess.run(
        [biov, "--config", str(config), "run", str(script)],
        env={
            **{k: v for k, v in os.environ.items() if not k.startswith("BIOV_")},
            # A stale inherited selection must not reach the script either.
            "BIOV_CONFIG": str(tmp_path / "missing.toml"),
        },
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == {
        "config": str(config.resolve()),
        "home": str(cache.resolve()),
    }
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize(
    "tool", ["unlisted-command", "conda:seqtk", "pypi:ruff", "npm:@scope/tool"]
)
@pytest.mark.parametrize("no_install", [False, True])
def test_exec_routes_any_command_through_configured_ssh(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, tool: str, no_install: bool
) -> None:
    """Keep remote paths and deployment details outside the invocation."""
    ssh_config = tmp_path / "ssh-config"
    ssh_config.write_text("Host compute\n  HostName example.invalid\n  Port 2222\n")
    selected = Settings.model_validate(
        {
            "execution_host": "compute",
            "ssh_config": ssh_config,
            "execution_cwd": Path("/remote/only/work"),
            "home": tmp_path / "must-not-be-created",
        }
    )
    monkeypatch.setattr(software, "settings", selected)
    calls = []

    def run(argv, **kwargs):
        calls.append(argv)
        assert kwargs == {"check": False}
        return subprocess.CompletedProcess(argv, 23)

    monkeypatch.setattr(subprocess, "run", run)
    arguments = [
        tool,
        "--help",
        "--cwd",
        "native-work",
        "--no-install",
        "a;$(touch forbidden)",
        "",
    ]
    options = ["--no-install"] if no_install else []
    result = CliRunner().invoke(app, ["exec", *options, *arguments])
    assert result.exit_code == 23, result.output
    argv = calls.pop()
    assert argv[:6] == ["ssh", "-T", "-F", str(ssh_config), "--", "compute"]
    remote = shlex.split(argv[-1])
    assert remote == [
        "biov",
        "exec",
        *options,
        "--cwd",
        "/remote/only/work",
        "--",
        *arguments,
    ]
    assert not any("BIOV_" in argument for argument in remote)
    assert not selected.home.exists()
    assert not calls
    assert "remote" not in {item.name for item in app.registered_commands}


def test_remote_preserves_shell_sensitive_arguments(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Exercise the internal SSH shell boundary locally with literal input."""
    native_run = subprocess.run
    results = []

    def run(argv, **kwargs):
        assert argv[:4] == ["ssh", "-T", "--", "test-host"]
        result = native_run(
            ["/bin/sh", "-c", argv[-1]], cwd=tmp_path, capture_output=True
        )
        results.append(result)
        return result

    monkeypatch.setattr(subprocess, "run", run)
    arguments = (
        "a b",
        "'quoted'",
        "$(touch forbidden)",
        "x;touch forbidden",
        "line\nbreak",
        "",
    )
    result = run_remote(
        "test-host",
        (
            sys.executable,
            "-c",
            "import json,sys; print(json.dumps(sys.argv[1:])); sys.exit(19)",
            *arguments,
        ),
    )
    assert result.returncode == 19
    assert json.loads(results[0].stdout) == list(arguments)
    assert not (tmp_path / "forbidden").exists()
    with pytest.raises(ValueError):
        run_remote("-oProxyCommand=bad", ("ls",))


def test_source_program_preserves_output_files_and_exit_status(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fake_uv
) -> None:
    """Run a source-delegated program, preserving output files and native status."""
    monkeypatch.setattr(software, "settings", Settings(home=tmp_path / "cache"))
    result = CliRunner().invoke(
        app,
        [
            "exec",
            "--cwd",
            str(tmp_path),
            "--",
            "pypi:python",
            "-c",
            "from pathlib import Path; import sys; Path('result.txt').write_text(sys.argv[1]); sys.exit(17)",
            "literal $value",
        ],
    )
    assert result.exit_code == 17, result.output
    assert (tmp_path / "result.txt").read_text() == "literal $value"
