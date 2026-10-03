"""Managed execution records, complete data reuse and bounded presentation."""

import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import unquote, urlsplit

import pytest

from biov import analysis


def local(reference: str) -> Path:
    """Return the test client's path for a local file URI."""
    return Path(unquote(urlsplit(reference).path))


@pytest.fixture
def runtime(monkeypatch, tmp_path):
    """Replace only deployment; execute the real worker and caller Python code.

    Returns:
        Mutable settings for the temporary local execution context.
    """
    config = SimpleNamespace(
        analysis_root=tmp_path / "results",
        analysis_base_url=None,
        execution_host=None,
    )
    monkeypatch.setattr(analysis, "settings", config)
    manifest = tmp_path / "pyproject.toml"
    manifest.write_text("[tool.pixi.environments]\npython = []\n")
    manifest.with_name("pixi.lock").write_text("version: 6\n")
    manager = tmp_path / "pixi"
    manager.write_text(f'#!{sys.executable}\nprint("pixi 0.81.0")\n')
    manager.chmod(0o700)
    monkeypatch.setattr(analysis, "pixi_command", lambda: (str(manager),))
    monkeypatch.setattr(analysis, "manifest_source", lambda: manifest)
    monkeypatch.setattr(analysis, "resolve_environment", lambda name: name)
    native_popen = subprocess.Popen

    def launch(command, **kwargs):
        if "from biov.cli import app; app()" in command:
            kwargs["stdout"].write(b"simulated Pixi preparation output\n")
            kwargs["stdout"].flush()
            return native_popen([sys.executable, *command[-2:]], **kwargs)
        return native_popen(command, **kwargs)

    monkeypatch.setattr(subprocess, "Popen", launch)
    return config


def test_complete_reuse_logs_and_failed_step_retention(runtime, capsys):
    """Reuse all rows and retain earlier success when a later check fails."""
    first = analysis.run_analysis(
        analysis.AnalysisRequest(
            name="source",
            environment="python",
            inputs={},
            outputs={"data.csv": "csv"},
            code="from pathlib import Path\nprint('program log')\nPath('data.csv').write_text('value\\n1\\n2\\n3\\n')\n",
        )
    )
    assert first["status"] == "succeeded"
    output = first["outputs"]["data.csv"]
    assert output["preview"]["rows"] == [[1], [2]]
    assert output["preview"]["total_rows"] is None
    assert output["preview"]["truncated"] is True
    assert (
        "program log"
        in analysis.read_analysis_file(first["logs"]["stdout"])[0].decode()
    )
    assert "program log" not in capsys.readouterr().out
    saved_first = local(first["record"]).read_bytes()

    second = analysis.run_analysis(
        analysis.AnalysisRequest(
            name="sum",
            environment="python",
            inputs={"values": output["uri"]},
            outputs={"sum.csv": "csv"},
            code="""import csv, json, sys
from pathlib import Path
inputs = json.loads(Path(sys.argv[1]).read_text())
with open(inputs['values']) as source:
    values = list(csv.DictReader(source))
Path('sum.csv').write_text('sum\\n' + str(sum(int(row['value']) for row in values)) + '\\n')
""",
        )
    )
    assert second["outputs"]["sum.csv"]["preview"]["rows"] == [[6]]
    failed = analysis.run_analysis(
        analysis.AnalysisRequest(
            name="failed",
            environment="python",
            inputs={"values": output["uri"]},
            outputs={"partial.csv": "csv"},
            code="""from pathlib import Path
import json
Path('partial.csv').write_text('partial\\n1\\n')
Path('checks.json').write_text(json.dumps([{'name': 'residues', 'passed': False, 'detail': 'unsupported residue'}]))
raise ValueError('original scientific diagnostic')
""",
        )
    )
    assert failed["status"] == "failed" and failed["stage"] == "data_checks"
    assert failed["exit_code"] != 0 and failed["outputs"] == {}
    assert failed["checks"][0]["passed"] is False
    assert (
        "original scientific diagnostic"
        in analysis.read_analysis_file(failed["logs"]["stderr"])[0].decode()
    )
    assert local(first["record"]).read_bytes() == saved_first
    assert analysis.inspect_analysis(failed["record"])["status"] == "failed"
    partial = local(failed["record"]).parent / "partial.csv"
    with pytest.raises(ValueError, match="not a completed"):
        analysis.read_analysis_file(partial.as_uri())


def test_changed_saved_output_is_not_read_or_reused(runtime):
    """Reject results whose saved content was replaced or deleted."""
    result = analysis.run_analysis(
        analysis.AnalysisRequest(
            name="file",
            environment="python",
            outputs={"result.txt": "file"},
            code="from pathlib import Path\nPath('result.txt').write_text('original')\n",
        )
    )
    uri = result["outputs"]["result.txt"]["uri"]
    local(uri).write_text("replaced")
    with pytest.raises(ValueError, match="identity mismatch"):
        analysis.inspect_analysis(result["record"])
    with pytest.raises(ValueError, match="identity mismatch"):
        analysis.read_analysis_file(uri)
    reused = analysis.run_analysis(
        analysis.AnalysisRequest(
            name="reuse",
            environment="python",
            inputs={"source": uri},
            outputs={"unused": "file"},
            code="raise AssertionError('must not execute')",
        )
    )
    assert reused["status"] == "failed" and reused["stage"] == "data_checks"
    assert "identity mismatch" in reused["diagnostic"]
    local(uri).unlink()
    with pytest.raises(FileNotFoundError):
        analysis.inspect_analysis(result["record"])


def test_snapshot_identity_and_copy_changes(runtime, tmp_path, monkeypatch):
    """Detect changes during copying and later changes to the execution input."""
    source = tmp_path / "source.txt"
    source.write_text("original")
    request = analysis.AnalysisRequest(
        name="mutate",
        environment="python",
        inputs={"source": str(source)},
        outputs={"unused": "file"},
        code="""from pathlib import Path
import json, sys
path = Path(json.loads(Path(sys.argv[1]).read_text())['source'])
path.chmod(0o600)
path.write_text('mutated')
""",
    )
    changed = analysis.run_analysis(request)
    assert (
        changed["status"] == "failed"
        and "Execution input identity mismatch" in changed["diagnostic"]
    )
    assert source.read_text() == "original"
    copy = analysis.shutil.copyfileobj

    def changing_copy(original, destination):
        copy(original, destination)
        source.write_text("changed during copy")

    monkeypatch.setattr(analysis.shutil, "copyfileobj", changing_copy)
    changed = analysis.run_analysis(request)
    assert "Input changed while copying" in changed["diagnostic"]
    assert "command" not in json.loads(local(changed["record"]).read_text())


@pytest.mark.parametrize("confirmed", [False, True])
def test_interruption_keeps_known_facts_without_finalization(
    runtime, monkeypatch, confirmed
):
    """Keep unconfirmed execution unknown and preserve saved evidence."""
    native_popen = subprocess.Popen

    def interrupted(command, **kwargs):
        if "from biov.cli import app; app()" not in command:
            return native_popen(command, **kwargs)
        if confirmed:
            (Path(command[-1]) / "execution.json").write_text(
                json.dumps({"status": "running", "pid": 12345})
            )
        raise KeyboardInterrupt

    monkeypatch.setattr(subprocess, "Popen", interrupted)
    with pytest.raises(KeyboardInterrupt):
        analysis.run_analysis(
            analysis.AnalysisRequest(
                name="interrupted",
                environment="python",
                outputs={"unused": "file"},
                code="pass",
            )
        )
    record = next(runtime.analysis_root.glob("*/record.json"))
    before = record.read_bytes()
    result = analysis.inspect_analysis(record.as_uri())
    assert result["status"] == "unknown" and result["exit_code"] is None
    assert record.read_bytes() == before
    if confirmed:
        assert result["environment"]["execution"]["status"] == "running"


@pytest.mark.parametrize("failure", ["confirmation_write", "evidence_read"])
def test_record_failure_after_launch_keeps_execution_unknown(
    runtime, monkeypatch, failure
):
    """Do not report scientific failure when confirmation storage fails."""
    save = analysis._save
    read = Path.read_text
    launch = subprocess.Popen
    processes = []
    injected = False

    def observe(command, **kwargs):
        process = launch(command, **kwargs)
        if "from biov.cli import app; app()" in command:
            processes.append(process)
        return process

    def fail_save(path, record):
        nonlocal injected
        if (
            not injected
            and failure == "confirmation_write"
            and "launcher_pid" in record
        ):
            injected = True
            raise OSError("confirmation storage failed")
        save(path, record)

    def fail_read(path, *args, **kwargs):
        nonlocal injected
        if (
            not injected
            and failure == "evidence_read"
            and path.name == "execution.json"
        ):
            injected = True
            raise OSError("execution evidence unavailable")
        return read(path, *args, **kwargs)

    monkeypatch.setattr(subprocess, "Popen", observe)
    monkeypatch.setattr(analysis, "_save", fail_save)
    monkeypatch.setattr(Path, "read_text", fail_read)
    result = analysis.run_analysis(
        analysis.AnalysisRequest(
            name="storage-failure",
            environment="python",
            outputs={"result.txt": "file"},
            code="from pathlib import Path\nPath('result.txt').write_text('finished')\n",
        )
    )
    for process in processes:
        process.wait(timeout=10)
    assert injected and result["status"] == "unknown"
    assert result["exit_code"] is None and result["outputs"] == {}
    inspected = analysis.inspect_analysis(result["record"])
    assert inspected["status"] == "unknown"
    assert inspected["environment"]["execution"]["status"] == "completed"


def test_source_can_change_after_snapshot_but_recorded_output_cannot(
    runtime, tmp_path, monkeypatch
):
    """Use copied bytes and reject a reused result changed before its copy."""
    source = tmp_path / "original.txt"
    source.write_text("snapshot bytes")
    result = analysis.run_analysis(
        analysis.AnalysisRequest(
            name="stable-copy",
            environment="python",
            inputs={"source": str(source)},
            outputs={"result.txt": "file"},
            parameters={"original": str(source)},
            code="""from pathlib import Path
import json, sys
inputs, parameters = [json.loads(Path(name).read_text()) for name in sys.argv[1:]]
Path(parameters['original']).write_text('source updated')
Path('result.txt').write_bytes(Path(inputs['source']).read_bytes())
""",
        )
    )
    assert result["status"] == "succeeded"
    uri = result["outputs"]["result.txt"]["uri"]
    assert analysis.read_analysis_file(uri)[0] == b"snapshot bytes"
    assert source.read_text() == "source updated"
    snapshot = analysis._snapshot

    def changed_before_copy(source, target):
        source.write_text("changed after reference verification")
        return snapshot(source, target)

    monkeypatch.setattr(analysis, "_snapshot", changed_before_copy)
    reused = analysis.run_analysis(
        analysis.AnalysisRequest(
            name="identity-race",
            environment="python",
            inputs={"source": uri},
            outputs={"result.txt": "file"},
            code="raise AssertionError('must not execute')",
        )
    )
    assert (
        reused["status"] == "failed"
        and "changed before snapshot" in reused["diagnostic"]
    )


def test_concurrent_runs_bounded_previews_and_complete_resource_limit(runtime):
    """Isolate runs, bound long previews and provide a full-file retrieval route."""
    request = analysis.AnalysisRequest(
        name="bounded",
        environment="python",
        outputs={"wide.csv": "csv", "sequence.fa": "fasta", "large.bin": "file"},
        code="""from pathlib import Path
Path('wide.csv').write_text('text\\n' + 'x' * 50000 + '\\n')
Path('sequence.fa').write_text('>id ' + 'description' * 5000 + '\\n' + 'A' * 1000 + '\\n')
Path('large.bin').write_bytes(b'x' * 1048577)
""",
    )
    with ThreadPoolExecutor(max_workers=2) as pool:
        first, second = list(pool.map(analysis.run_analysis, [request, request]))
    assert first["status"] == second["status"] == "succeeded"
    assert first["record"] != second["record"]
    assert (
        len(json.dumps(first, ensure_ascii=False).encode())
        < analysis.MAX_RESPONSE_BYTES
    )
    assert first["outputs"]["wide.csv"]["preview"]["truncated"]
    assert first["outputs"]["sequence.fa"]["preview"]["truncated"]
    uri = first["outputs"]["large.bin"]["uri"]
    runtime.analysis_base_url = "http://127.0.0.1:8000/files"
    with pytest.raises(ValueError, match=r"download from http://127\.0\.0\.1"):
        analysis.read_analysis_file(uri)
    assert analysis.inspect_analysis(first["record"])["outputs"]["large.bin"][
        "download_url"
    ].startswith(runtime.analysis_base_url)
    saved = json.loads(local(first["record"]).read_text())
    assert saved["environment"]["manifest_sha256"]
    assert saved["environment"]["execution"]["interpreter"] == sys.executable
    assert saved["code"]["sha256"] == hashlib.sha256(request.code.encode()).hexdigest()


def test_large_conditions_and_checks_remain_in_record_with_explicit_summary_omissions(
    runtime,
):
    """Preserve complete conditions and checks while bounding their summary."""
    result = analysis.run_analysis(
        analysis.AnalysisRequest(
            name="conditions",
            environment="python",
            outputs={"result.txt": "file"},
            parameters={"long": "x" * 40000},
            requirements=["y" * 40000],
            code="""import json
from pathlib import Path
Path('result.txt').write_text('ok')
Path('checks.json').write_text(json.dumps([{'name': 'check', 'passed': True, 'detail': 'z' * 40000}]))
""",
        )
    )
    assert result["status"] == "succeeded"
    assert isinstance(result["checks"], list) and isinstance(
        result["requirements"], list
    )
    assert (
        len(json.dumps(result, ensure_ascii=False).encode())
        <= analysis.MAX_RESPONSE_BYTES - 512
    )
    assert {"checks", "requirements", "parameters"} <= set(result["omitted"])
    saved = json.loads(local(result["record"]).read_text())
    assert len(saved["parameters"]["long"]) == 40000
    assert len(saved["checks"][0]["detail"]) == 40000


@pytest.mark.parametrize(
    "name",
    ["../escape", "nested/file", "record.json", "STDOUT.LOG", "execution.json", "a..b"],
)
def test_request_rejects_reserved_or_ambiguous_output_names(name):
    """Reject names that could escape the run or overwrite engine records."""
    with pytest.raises(ValueError):
        analysis.AnalysisRequest(
            name="test", environment="python", code="pass", outputs={name: "file"}
        )
