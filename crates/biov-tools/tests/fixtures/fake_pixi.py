#!/usr/bin/env python3
"""Protocol fixture for native Rust tests; not a package manager."""

# ruff: noqa: T201, S101 - fixture deliberately exercises native stdio and assertions
import json
import pathlib
import sys

home = pathlib.Path(__file__).parent
args = sys.argv[1:]
with (home / "argv.jsonl").open("a") as out:
    out.write(json.dumps(args) + "\n")
if args == ["--version"]:
    print(
        (home / "version").read_text().strip()
        if (home / "version").exists()
        else "pixi 0.81.0"
    )
    sys.exit(0)
manifest = pathlib.Path(args[args.index("--manifest-path") + 1])
if args[0] == "info":
    print(
        json.dumps(
            {
                "platform": "linux-64",
                "version": "0.81.0",
                "environments_info": [
                    {
                        "name": name,
                        "prefix": str(
                            home / "external"
                            if (home / "redirect").exists()
                            else manifest.parent / ".pixi" / "envs" / name
                        ),
                    }
                    for name in ("samtools", "goatools")
                ],
            }
        )
    )
    sys.exit(0)
name = args[args.index("--environment") + 1]
prefix = manifest.parent / ".pixi" / "envs" / name
if args[0] == "install":
    if (home / "fail-install").exists():
        sys.exit(23)
    assert "--locked" in args
    (prefix / "conda-meta").mkdir(parents=True, exist_ok=True)
    (prefix / "bin").mkdir(exist_ok=True)
    (prefix / "bin" / name).write_text("#!/bin/sh\nexit 0\n")
    (prefix / "bin" / name).chmod(0o755)
    (prefix / "conda-meta" / "pixi").write_text(
        json.dumps(
            {
                "environment_name": name,
                "pixi_version": "0.81.0",
                "manifest_path": str(manifest),
                "resolved_platform": {"subdir": "linux-64"},
                "environment_lock_file_hash": "fixture",
            }
        )
    )
    sys.exit(0)
if args[0] == "run":
    assert "--as-is" in args
    native = args[args.index("--") + 1 :]
    (home / "native.json").write_text(
        json.dumps({"argv": native, "cwd": str(pathlib.Path.cwd())})
    )
    print("native stdout")
    print("native stderr", file=sys.stderr)
    sys.exit(37)
sys.exit(99)
