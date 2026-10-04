"""Verify the unified wheel's mandatory native executable and Python payload."""

# This build-verification script deliberately asserts distribution invariants.
# ruff: noqa: S101

import base64
import configparser
import csv
import email
import hashlib
import io
import sys
import zipfile
from pathlib import Path


def main(filename: str) -> None:
    """Require one binary, one ABI3 extension, complete assets and valid RECORD."""
    with zipfile.ZipFile(filename) as wheel:
        names = wheel.namelist()
        assert len(names) == len(set(names)), "wheel has duplicate archive members"
        binaries = [name for name in names if ".data/scripts/" in name]
        assert len(binaries) == 1, binaries
        binary = binaries[0]
        assert binary.endswith(("/biov", "/biov.exe")), binary
        executable = wheel.read(binary)
        assert executable.startswith(
            (
                b"\x7fELF",
                b"MZ",
                b"\xcf\xfa\xed\xfe",
                b"\xfe\xed\xfa\xcf",
                b"\xca\xfe\xba\xbe",
            )
        ), "biov must be a native executable"
        if not binary.endswith(".exe"):
            assert (wheel.getinfo(binary).external_attr >> 16) & 0o111, (
                "biov is not executable"
            )
        extensions = [
            name
            for name in names
            if name.startswith("biov/_native.") and name.endswith((".so", ".pyd"))
        ]
        assert len(extensions) == 1, extensions
        assert ".abi3." in extensions[0] or extensions[0].endswith(".pyd"), extensions
        entries = [
            name for name in names if name.endswith(".dist-info/entry_points.txt")
        ]
        assert len(entries) == 1, entries
        config = configparser.ConfigParser()
        config.read_string(wheel.read(entries[0]).decode())
        assert "biov" not in config["console_scripts"], "Python must not own biov"
        for package in ("biov", "crisprprimer"):
            for source in Path("src", package).rglob("*"):
                if not source.is_file() or "__pycache__" in source.parts:
                    continue
                if source.suffix in {".pyc", ".so", ".pyd", ".dylib"}:
                    continue
                relative = source.relative_to("src").as_posix()
                assert wheel.read(relative) == source.read_bytes(), relative
        records = [name for name in names if name.endswith(".dist-info/RECORD")]
        assert len(records) == 1, records
        rows = list(csv.reader(io.StringIO(wheel.read(records[0]).decode())))
        assert {row[0] for row in rows} == set(names), "RECORD member set differs"
        assert len(rows) == len(names), "RECORD has duplicate members"
        for name, digest, size in rows:
            if name == records[0]:
                assert digest == size == ""
                continue
            payload = wheel.read(name)
            expected = (
                base64.urlsafe_b64encode(hashlib.sha256(payload).digest())
                .rstrip(b"=")
                .decode()
            )
            assert digest == f"sha256={expected}", name
            assert size == str(len(payload)), name
        wheel_metadata = [name for name in names if name.endswith(".dist-info/WHEEL")]
        assert len(wheel_metadata) == 1, wheel_metadata
        metadata = email.message_from_bytes(wheel.read(wheel_metadata[0]))
        assert metadata["Root-Is-Purelib"] == "false"
        assert all(tag.startswith("cp312-abi3-") for tag in metadata.get_all("Tag", []))
        assert metadata.get_all("Tag"), "wheel has no compatibility tag"
        assert "-cp312-abi3-" in Path(filename).name, filename


if __name__ == "__main__":
    main(sys.argv[1])
