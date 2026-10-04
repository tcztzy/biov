#!/usr/bin/env python3
"""Protocol fixture for native Rust tests; not a package manager."""

# ruff: noqa: T201, S101 - fixture deliberately exercises native stdio and assertions
import json
import os
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
if args[0] in {"install", "reinstall"}:
    if (home / "fail-install").exists() or (
        args[0] == "reinstall" and (home / "fail-reinstall").exists()
    ):
        sys.exit(23)
    assert "--locked" in args
    if args[0] == "install" and (prefix / "conda-meta" / "pixi").exists():
        sys.exit(0)
    (prefix / "conda-meta").mkdir(parents=True, exist_ok=True)
    (prefix / "bin").mkdir(exist_ok=True)
    native_program = f"""#!/usr/bin/env python3
import json, pathlib, sys
home = pathlib.Path({str(home)!r})
(home / 'native.json').write_text(json.dumps({{'argv': [str(pathlib.Path(__file__)), *sys.argv[1:]], 'cwd': str(pathlib.Path.cwd())}}))
print('native stdout')
print('native stderr', file=sys.stderr)
sys.exit(37)
"""
    (prefix / "bin" / name).write_text(native_program)
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
if args[0] == "shell-hook":
    assert "--as-is" in args and "--json" in args
    activation_prefix = (
        home / "outside" if (home / "activation-redirect").exists() else prefix
    )
    print(
        json.dumps(
            {
                "environment_variables": {
                    "CONDA_PREFIX": str(activation_prefix),
                    "PIXI_PROJECT_MANIFEST": str(manifest),
                    "PIXI_PROJECT_ROOT": str(manifest.parent),
                    "PIXI_ENVIRONMENT_NAME": name,
                    "PATH": str(prefix / "bin") + os.pathsep + os.environ["PATH"],
                },
                "activation_scripts": [],
            }
        )
    )
    sys.exit(0)
sys.exit(99)
