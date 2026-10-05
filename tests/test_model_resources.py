"""Installed CLI and opt-in public selected-model-file portability acceptance.

BIOV_TEST_BINARY selects an independently installed native BioV executable.
BIOV_TEST_REAL_MODELS=1 enables the 869-byte public configuration download via
official hf/uv. Offline checks run under the existing child-only network filter;
the independent consumer uses isolated standard-library Python without BioV.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BINARY = os.environ.get("BIOV_TEST_BINARY", "")
REVISION = "f171d7baecaf37b5da5a3616d8833b9969753535"
FILES = ["config.json", "tokenizer_config.json"]
pytestmark = pytest.mark.skipif(
    not BINARY or sys.platform != "linux",
    reason="requires a separately installed Linux native BioV tool",
)


def _run(arguments, cwd, environment, *, offline=False):
    """Return a bounded native or independent-reader acceptance result."""
    if offline:
        arguments = [
            sys.executable,
            "-I",
            str(ROOT / "tests/_native_storage_sandbox.py"),
            *arguments,
        ]
    return subprocess.run(
        arguments,
        cwd=cwd,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )


def _record(result):
    """Require a successful command with exactly one JSON stdout result.

    Returns:
        The complete decoded command result.
    """
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def _verify_moved_bundle(bundle, cwd):
    """Verify the moved closure with no managers, network or BioV import.

    Returns:
        The moved path, inspected resource and independent-reader summary.
    """
    moved = cwd / "moved selected files"
    bundle.rename(moved)
    environment = {"PATH": "", "PYTHONPATH": "", "PYTHONNOUSERSITE": "1"}
    inspected = _record(
        _run([BINARY, "model", "inspect", str(moved)], cwd, environment, offline=True)
    )
    assert inspected["status"] == "verified"
    readme = (moved / "BIOV_MODEL_RESOURCE_README.md").read_text()
    script = readme.split("<<'PY'\n", 1)[1].split("\nPY\n", 1)[0]
    script += (
        "\nimport importlib.util\n"
        "require(importlib.util.find_spec('biov') is None, 'BioV unexpectedly available')\n"
    )
    reader = _run(
        [sys.executable, "-I", "-S", "-c", script, str(moved)],
        cwd,
        environment,
        offline=True,
    )
    summary = _record(reader)
    assert summary["verified_files"] == len(FILES)
    return moved, inspected, summary


def test_installed_selected_model_files_and_offline_reader(tmp_path):
    """Installed packaging retains native routing and portable complete records."""
    hf = tmp_path / "hf"
    shutil.copyfile(ROOT / "crates/biov-cli/tests/fixtures/fake_hf.py", hf)
    hf.chmod(0o755)
    bundle = tmp_path / "selected files"
    environment = {**os.environ, "BIOV_FAKE_CLIENT_HOME": str(tmp_path)}
    for name in ("HF_ENDPOINT", "HUGGINGFACE_CO_STAGING", "BIOV_HF_BIN"):
        environment.pop(name, None)
    downloaded = _record(
        _run(
            [
                BINARY,
                "model",
                "download",
                "--revision",
                REVISION,
                "--local-dir",
                str(bundle),
                "--hf",
                str(hf),
                "example/model",
                *FILES,
            ],
            tmp_path,
            environment,
        )
    )
    assert downloaded["status"] == "downloaded"
    assert downloaded["resource"]["scope"] == "selected_files"
    moved, inspected, summary = _verify_moved_bundle(bundle, tmp_path)
    assert inspected["resource"] == downloaded["resource"]
    assert summary["complete_payload_bytes"] == sum(
        item["bytes"] for item in downloaded["resource"]["inventory"]
    )
    (moved / "config.json").write_bytes(b"damaged")
    rejected = _run(
        [BINARY, "model", "inspect", str(moved)],
        tmp_path,
        {"PATH": ""},
        offline=True,
    )
    assert rejected.returncode == 2
    assert rejected.stdout == ""


@pytest.mark.skipif(
    os.environ.get("BIOV_TEST_REAL_MODELS") != "1",
    reason="real public model-resource download is explicitly opt-in",
)
def test_real_public_configuration_download_reuse_and_moved_reader(tmp_path):
    """Official acquisition yields complete native JSON without model execution."""
    uv = os.environ.get("BIOV_TEST_REAL_UV") or shutil.which("uv")
    assert uv, "real model acceptance requires the official uv executable"
    bundle = tmp_path / "public configuration"
    environment = {**os.environ, "HF_HUB_DISABLE_IMPLICIT_TOKEN": "1"}
    for name in ("HF_TOKEN", "HF_ENDPOINT", "HUGGINGFACE_CO_STAGING", "BIOV_HF_BIN"):
        environment.pop(name, None)
    for name in ("HF_HOME", "UV_TOOL_DIR", "UV_TOOL_BIN_DIR", "UV_CACHE_DIR"):
        environment[name] = str(tmp_path / name.lower())
    arguments = [
        BINARY,
        "model",
        "download",
        "--revision",
        REVISION,
        "--local-dir",
        str(bundle),
        "--uv",
        str(uv),
        "hf-internal-testing/tiny-random-bert",
        *FILES,
    ]
    downloaded = _record(_run(arguments, tmp_path, environment))
    assert downloaded["status"] == "downloaded"
    inventory = downloaded["resource"]["inventory"]
    assert [(item["path"], item["bytes"]) for item in inventory] == [
        ("config.json", 548),
        ("tokenizer_config.json", 321),
    ]
    reused = _record(_run(arguments, tmp_path, {"PATH": ""}, offline=True))
    assert reused["status"] == "reused"
    assert reused["resource"] == downloaded["resource"]
    moved, inspected, summary = _verify_moved_bundle(bundle, tmp_path)
    assert inspected["resource"] == downloaded["resource"]
    assert summary["complete_payload_bytes"] == 869
    # Read every native JSON field with an ordinary reader, without model imports.
    native_reader = (
        "import importlib.util,json,sys; from pathlib import Path; "
        "assert importlib.util.find_spec('biov') is None; "
        "root=Path(sys.argv[1]); "
        "config=json.loads((root/'config.json').read_bytes()); "
        "tokenizer=json.loads((root/'tokenizer_config.json').read_bytes()); "
        "print(json.dumps({'model_type':config['model_type'],"
        "'configuration_fields':len(config),'tokenizer_fields':len(tokenizer)}))"
    )
    native = _record(
        _run(
            [sys.executable, "-I", "-S", "-c", native_reader, str(moved)],
            tmp_path,
            {"PATH": ""},
            offline=True,
        )
    )
    assert native["model_type"] == "bert"
    assert native["configuration_fields"] > 0 and native["tokenizer_fields"] > 0
