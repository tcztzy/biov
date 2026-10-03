"""Check native source distribution boundaries and required build inputs."""


# This build-verification script deliberately asserts distribution invariants.
# ruff: noqa: S101

import re
import sys
import tarfile
from pathlib import Path, PurePosixPath


def main(filename: str) -> None:
    """Reject missing sources, build artifacts and workstation-specific paths."""
    with tarfile.open(filename) as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        paths = {member.name.split("/", 1)[1] for member in members}
        required = {
            "Cargo.toml",
            "Cargo.lock",
            "uv.lock",
            "rust-toolchain.toml",
            "pyproject.toml",
            "LICENSE",
            "README.md",
            "SPEC.md",
            "src/biov/_native.pyi",
            "src/biov/py.typed",
            "crates/biov-core/src/sequence.rs",
            "crates/biov-python/src/lib.rs",
            "crates/biov-core/src/bin/biov-core.rs",
            "tests/test_native_sequence.py",
            "docs/guides/sequence-contract.md",
        }
        assert required <= paths, required - paths
        forbidden = {"target", ".venv", "__pycache__", ".git", ".pixi", "site", "dist"}
        selected = {
            "AGENTS.md",
            "Cargo.toml",
            "Cargo.lock",
            "uv.lock",
            "LICENSE",
            "README.md",
            "SPEC.md",
            "mkdocs.yml",
            "pyproject.toml",
            "rust-toolchain.toml",
        }
        for directory in (
            "src/biov",
            "src/crisprprimer",
            "crates",
            "docs",
            "tests",
            "scripts",
        ):
            selected.update(
                path.as_posix()
                for path in Path(directory).rglob("*")
                if path.is_file()
                and not forbidden.intersection(path.parts)
                and path.suffix not in {".so", ".pyc", ".whl"}
            )
        selected.discard("tests/test_plugin_distribution.py")
        assert paths - {"PKG-INFO"} == selected, (
            f"sdist mismatch: missing={sorted(selected - paths)}, "
            f"unexpected={sorted(paths - selected - {'PKG-INFO'})}"
        )
        workstation = re.compile(rb"/(?:Users|Volumes)/")
        for member in members:
            path = PurePosixPath(member.name)
            assert not forbidden.intersection(path.parts), path
            assert path.suffix not in {".so", ".pyc", ".whl"}, path
            stream = archive.extractfile(member)
            assert stream is not None
            assert not workstation.search(stream.read()), path


if __name__ == "__main__":
    main(sys.argv[1])
