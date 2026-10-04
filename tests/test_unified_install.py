"""Acceptance against a real wheel/sdist installed by uv tool outside checkout.

Set BIOV_TEST_BINARY to the public native tool path to run these checks. The
normal source suite skips them, since cargo artifacts cannot establish Python
packaging or paired-interpreter behavior.

BIOV_TEST_REAL_TOOLS=1 also requires BIOV_TEST_WHEEL and BIOV_TEST_REAL_PIXI.
That packaging gate uses actual Pixi/uv installs, runs native scientific checks,
and removes a separate management installation before repeating those checks.
"""

import json
import os
import selectors
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

BINARY = os.environ.get("BIOV_TEST_BINARY", "")
pytestmark = pytest.mark.skipif(
    not BINARY, reason="requires a separately installed native BioV tool"
)


def test_native_help_and_legacy_help_without_ambient_python(tmp_path):
    """Installed native routes and the explicit legacy bridge ignore PATH Python."""
    fake = tmp_path / "biov"
    fake.mkdir()
    (fake / "__init__.py").write_text("raise RuntimeError('stale checkout imported')")
    environment = {**os.environ, "PATH": "", "PYTHONPATH": str(tmp_path)}
    native = subprocess.run(
        [BINARY, "--help"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert native.returncode == 0
    assert native.stdout == ""
    assert "mcp-native" in native.stderr
    assert "biov-rs" not in native.stderr
    python = subprocess.run(
        [BINARY, "python", "--help"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert python.returncode == 0, python.stderr
    assert "inspect-analysis" in python.stdout
    assert "setup" in python.stdout


def test_python_run_uses_paired_interpreter_and_preserves_argv_status(tmp_path):
    """A real installed CLI launches a normal script in its own tool environment."""
    script = tmp_path / "script.py"
    script.write_text(
        "import json, sys, biov\nprint(json.dumps({'interpreter':sys.executable,'package':biov.__file__,'args':sys.argv[1:]}))\nraise SystemExit(37)\n"
    )
    arguments = ["", "O'Connor", "$(touch must-not-exist)", "--native"]
    result = subprocess.run(
        [BINARY, "python", "run", str(script), "--", *arguments],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 37, result.stderr
    record = json.loads(result.stdout)
    executable = Path(BINARY).resolve()
    assert Path(record["interpreter"]).parent == executable.parent
    assert record["args"] == arguments
    assert "site-packages" in record["package"]
    assert not (tmp_path / "must-not-exist").exists()


def _read_rpc(process: subprocess.Popen[str]) -> dict:
    """Read an actual protocol response with a bounded wait.

    Returns:
        The complete JSON-RPC response.
    """
    assert process.stdout is not None
    selector = selectors.DefaultSelector()
    try:
        selector.register(process.stdout, selectors.EVENT_READ)
        assert selector.select(timeout=30), "MCP produced no response within 30 seconds"
        line = process.stdout.readline()
        assert line, "MCP ended before returning a response"
        return json.loads(line)
    finally:
        selector.close()


@pytest.mark.skipif(
    sys.platform == "win32", reason="pipe selector acceptance is POSIX-specific"
)
@pytest.mark.parametrize(
    "route,server_name,required_tool",
    [
        ("mcp", "biov", "parse_identifiers"),
        ("mcp-native", "biov-native", "dataset_open"),
    ],
)
def test_distinct_legacy_and_native_mcp_stdio(
    tmp_path, route, server_name, required_tool
):
    """Bare mcp retains Python tools; explicit mcp-native retains Rust datasets."""
    args = [BINARY, route]
    if route == "mcp-native":
        args += ["--data-root", str(tmp_path), "--output-root", str(tmp_path)]
    process = subprocess.Popen(
        args,
        cwd=tmp_path,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert process.stdin is not None and process.stdout is not None
    try:
        initialize = {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": "2025-06-18",
                "capabilities": {},
                "clientInfo": {"name": "unified-install-test", "version": "1"},
            },
        }
        process.stdin.write(json.dumps(initialize) + "\n")
        process.stdin.flush()
        result = _read_rpc(process)
        assert result["id"] == 1
        assert result["result"]["serverInfo"]["name"] == server_name
        process.stdin.write(
            json.dumps({"jsonrpc": "2.0", "method": "notifications/initialized"}) + "\n"
        )
        process.stdin.write(
            json.dumps(
                {"jsonrpc": "2.0", "id": 2, "method": "tools/list", "params": {}}
            )
            + "\n"
        )
        process.stdin.flush()
        response = _read_rpc(process)
        assert response["id"] == 2
        names = {tool["name"] for tool in response["result"]["tools"]}
        assert required_tool in names
        assert ("dataset_open" in names) == (route == "mcp-native")
        process.stdin.close()
        assert process.wait(timeout=30) == 0
        assert process.stdout.read() == ""
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=10)


def _checked_command(command, *, cwd, environment, evidence, label, timeout=600):
    """Run real upstream commands and retain complete process evidence.

    Returns:
        The successful completed process, including unabridged stdout/stderr.
    """
    result = subprocess.run(
        [str(arg) for arg in command],
        cwd=cwd,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout,
    )
    evidence.mkdir(exist_ok=True)
    (evidence / f"{label}.stdout").write_text(result.stdout)
    (evidence / f"{label}.stderr").write_text(result.stderr)
    (evidence / f"{label}.json").write_text(
        json.dumps(
            {"argv": [str(arg) for arg in command], "exitCode": result.returncode},
            indent=2,
        )
        + "\n"
    )
    assert result.returncode == 0, (
        f"{label}: {command!r}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    return result


@pytest.mark.skipif(
    os.environ.get("BIOV_TEST_REAL_TOOLS") != "1",
    reason="explicit real Pixi/uv lifecycle gate requires package cache or network",
)
def test_real_upstream_tools_survive_uninstalling_management_wheel(tmp_path):
    """Use real managers, then run scientific commands with BioV uninstalled.

    This opt-in gate is mandatory in packaging CI for both the checkout wheel
    and the wheel rebuilt from its unpacked source distribution. No manager,
    executable, package inventory, or scientific output is mocked here.
    """
    wheel_setting = os.environ.get("BIOV_TEST_WHEEL")
    pixi_setting = os.environ.get("BIOV_TEST_REAL_PIXI")
    uv_setting = os.environ.get("BIOV_TEST_REAL_UV") or shutil.which("uv")
    assert wheel_setting, "BIOV_TEST_WHEEL is required by the enabled real gate"
    assert pixi_setting, "BIOV_TEST_REAL_PIXI is required by the enabled real gate"
    assert uv_setting, "uv is required by the enabled real gate"
    wheel = Path(wheel_setting).resolve(strict=True)
    pixi = Path(pixi_setting).resolve(strict=True)
    uv = Path(uv_setting).resolve(strict=True)
    checkout = Path(__file__).resolve().parents[1]

    # The proof owns another real management install. Removing it cannot remove
    # BIOV_TEST_BINARY or introduce test-order dependencies in the other gates.
    cwd = tmp_path / "unrelated"
    cwd.mkdir()
    copied_wheel = cwd / wheel.name
    shutil.copyfile(wheel, copied_wheel)
    management_tools = tmp_path / "management-tools"
    management_bin = tmp_path / "management-bin"
    root = tmp_path / "upstream tools 'quoted'"
    evidence = cwd / "evidence"
    environment = {
        **os.environ,
        "UV_TOOL_DIR": str(management_tools),
        "UV_TOOL_BIN_DIR": str(management_bin),
        "PYTHONPATH": "",
        "PIXI_COLOR": "never",
    }

    def run(command, label, *, env=None):
        return _checked_command(
            command,
            cwd=cwd,
            environment=environment if env is None else env,
            evidence=evidence,
            label=label,
        )

    run([pixi, "--version"], "pixi-version")
    run([uv, "--version"], "uv-version")
    run(
        [
            uv,
            "tool",
            "install",
            "--python",
            os.environ.get("BIOV_TEST_TOOL_PYTHON", "3.12"),
            copied_wheel,
        ],
        "install-management-wheel",
    )
    manager = management_bin / "biov"
    assert manager.is_file()
    assert not (management_bin / "biov-rs").exists()
    management_prefix = manager.resolve().parent.parent
    assert manager.resolve() != Path(BINARY).resolve()
    paired_base = run(
        [
            management_prefix / "bin" / "python",
            "-I",
            "-c",
            "import pathlib,sys; print(pathlib.Path(sys._base_executable).resolve())",
        ],
        "management-base-interpreter",
    ).stdout.strip()
    options = ["--environment-root", root, "--pixi", pixi, "--uv", uv]
    empty = run([manager, "list", *options], "empty-inventory")
    assert "No native tools installed" in empty.stdout
    assert not root.exists(), "empty inventory must not provision environments"
    for name in ("samtools", "goatools"):
        installed = run([manager, "install", *options, name], f"install-{name}")
        assert "export PATH=" in installed.stderr

    pixi_bin = root / "pixi-global" / "bin"
    uv_bin = root / "uv-tools" / "bin"
    upstream_tools = root / "uv-tools" / "tools"
    tool_environment = {
        **environment,
        "PATH": os.pathsep.join([str(pixi_bin), str(uv_bin), environment["PATH"]]),
    }
    for command, directory in (
        ("samtools", pixi_bin),
        ("goatools", uv_bin),
        ("find_enrichment.py", uv_bin),
    ):
        assert shutil.which(command, path=tool_environment["PATH"]) == str(
            directory / command
        )
        assert os.access(directory / command, os.X_OK)
    assert not (root / "runners").exists()
    assert not (root / "installed").exists()
    assert not list(root.rglob("biov")), "upstream tools must not retain a BioV runner"

    listing = run([manager, "list", *options], "installed-inventory")
    assert "Pixi global (" in listing.stdout
    assert "uv tool (" in listing.stdout
    assert "samtools" in listing.stdout and "1.24" in listing.stdout
    assert "goatools" in listing.stdout and "1.6.5" in listing.stdout
    goat_python = upstream_tools / "goatools" / "bin" / "python"
    assert goat_python.resolve() == Path(paired_base)
    assert not goat_python.resolve().is_relative_to(management_prefix)
    packages = run(
        [
            goat_python,
            "-I",
            "-c",
            (
                "import importlib.metadata as m, importlib.util, json; "
                "print(json.dumps({p:m.version(p) for p in ['goatools','statsmodels']})); "
                "assert importlib.util.find_spec('biov') is None"
            ),
        ],
        "goatools-package-versions",
    )
    assert json.loads(packages.stdout) == {
        "goatools": "1.6.5",
        "statsmodels": "0.14.6",
    }

    # Copy the standalone acceptance script and native inputs; no analysis or
    # tool command below depends on source-checkout imports or relative paths.
    (cwd / "scripts").mkdir()
    script = cwd / "scripts" / "validate_goatools.py"
    shutil.copyfile(checkout / "scripts" / script.name, script)
    shutil.copytree(
        checkout / "tests" / "fixtures" / "goatools",
        cwd / "tests" / "fixtures" / "goatools",
    )
    sam = cwd / "reads.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        "@SQ\tSN:chr1\tLN:100\n"
        "read1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        "read2\t4\t*\t0\t0\t*\t*\t0\t0\tTGCA\tIIII\n"
    )
    data = root / "biological-data"
    data.mkdir()
    retained = data / "retain.txt"
    retained.write_text("Caller-owned scientific data\n")
    locked_marker = root / "workspaces" / "existing-locked-workflow" / "keep.txt"
    locked_marker.parent.mkdir(parents=True)
    locked_marker.write_text("Existing separate locked-workflow content\n")

    def direct_science(label):
        versions = run(
            ["samtools", "--version"], f"{label}-samtools-version", env=tool_environment
        )
        assert "samtools 1.24" in versions.stdout
        run(["goatools", "--help"], f"{label}-goatools-help", env=tool_environment)
        run(
            ["find_enrichment.py", "--help"],
            f"{label}-enrichment-help",
            env=tool_environment,
        )
        bam = cwd / f"{label}.bam"
        run(
            ["samtools", "view", "-b", "-o", bam, sam],
            f"{label}-sam-to-bam",
            env=tool_environment,
        )
        count = run(
            ["samtools", "view", "-c", bam],
            f"{label}-alignment-count",
            env=tool_environment,
        )
        assert count.stdout.strip() == "2"
        flags = run(
            ["samtools", "flagstat", bam], f"{label}-flagstat", env=tool_environment
        )
        assert "2 + 0 in total" in flags.stdout
        assert "1 + 0 mapped" in flags.stdout
        output = cwd / f"{label}-goatools"
        run(
            [sys.executable, "-I", script, output],
            f"{label}-go-enrichment",
            env=tool_environment,
        )
        verification = json.loads((output / "validation.json").read_text())
        assert verification["verified"] is True
        assert verification["exactExpectedLeafPvalues"] == {
            "uncorrected": "1/210",
            "bonferroni": "1/70",
            "fdr_bh": "1/140",
        }
        return output

    before = direct_science("before-removal")
    # Minimal uninstall now follows each backend's native environment removal.
    # User input/output and the other backend's installation stay untouched.
    sam_prefix = root / "pixi-global" / "envs" / "samtools"
    assert sam_prefix.is_dir()
    run([manager, "uninstall", *options, "samtools"], "uninstall-samtools")
    assert not (pixi_bin / "samtools").exists()
    assert not sam_prefix.exists()
    assert (uv_bin / "goatools").is_file()
    run([manager, "install", *options, "samtools"], "reinstall-samtools")
    run([manager, "uninstall", *options, "goatools"], "uninstall-goatools")
    assert not (uv_bin / "goatools").exists()
    assert not (upstream_tools / "goatools").exists()
    assert (pixi_bin / "samtools").is_file()
    run([manager, "install", *options, "goatools"], "reinstall-goatools")
    assert retained.read_text() == "Caller-owned scientific data\n"
    assert locked_marker.read_text() == "Existing separate locked-workflow content\n"
    assert (before / "results.tsv").is_file()

    run([uv, "tool", "uninstall", "biov"], "uninstall-management-wheel")
    assert not manager.exists()
    assert not management_prefix.exists()
    assert not (root / "runners").exists()
    after = direct_science("after-removal")
    assert (before / "results.tsv").read_bytes() == (after / "results.tsv").read_bytes()
    assert retained.read_text() == "Caller-owned scientific data\n"
    assert locked_marker.read_text() == "Existing separate locked-workflow content\n"
