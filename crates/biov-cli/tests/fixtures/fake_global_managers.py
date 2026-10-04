#!/usr/bin/env python3
"""CLI protocol fixture only; never evidence of installation or scientific use."""

# ruff: noqa: T201, S101 - subprocess fixture records arguments and writes test output
import json
import os
import shutil
import sys
from pathlib import Path

arguments = sys.argv[1:]
with Path(os.environ["BIOV_TEST_GLOBAL_LOG"]).open("a") as stream:
    stream.write(
        json.dumps(
            {
                "argv": arguments,
                **{
                    key: os.environ.get(key)
                    for key in ("PIXI_HOME", "UV_TOOL_DIR", "UV_TOOL_BIN_DIR")
                },
            }
        )
        + "\n"
    )
if arguments == ["--version"]:
    print("pixi 0.81.0" if "pixi" in Path(sys.argv[0]).name else "uv 0.12.19")
    sys.exit(0)
assert arguments[0] in ("global", "tool"), arguments
native = arguments[0] == "global"
root = Path(os.environ["PIXI_HOME"] if native else os.environ["UV_TOOL_DIR"])
bin_directory = root / "bin" if native else Path(os.environ["UV_TOOL_BIN_DIR"])
name = "samtools" if native else "goatools"
prefix = root / "envs" / name if native else root / name
action = arguments[1]
if action == "install":
    prefix.mkdir(parents=True, exist_ok=True)
    bin_directory.mkdir(parents=True, exist_ok=True)
    executable = bin_directory / name
    payload = "#!/bin/sh\nprintf 'protocol fixture only\\n'\n"
    if native:
        trampoline = bin_directory / "trampoline_configuration" / "trampoline_bin"
        trampoline.parent.mkdir(exist_ok=True)
        trampoline.write_text(payload)
        executable.write_text(payload)
        executable.chmod(0o755)
    else:
        (prefix / "bin").mkdir()
        entry = prefix / "bin" / name
        entry.write_text(payload)
        entry.chmod(0o755)
        executable.symlink_to(entry)
elif action == "list":
    if native and "--json" in arguments:
        print(
            json.dumps(
                [
                    {
                        "name": "samtools",
                        "dependencies": [{"name": "samtools", "version": "1.24"}],
                        "exposed": [
                            {"exposed_name": "samtools", "executable": "samtools"}
                        ],
                    }
                ]
                if prefix.exists()
                else []
            )
        )
    elif prefix.exists():
        print("samtools 1.24" if native else "goatools v1.6.5\n- goatools")
elif action == "uninstall":
    shutil.rmtree(prefix)
    (bin_directory / name).unlink()
else:
    raise AssertionError(arguments)
