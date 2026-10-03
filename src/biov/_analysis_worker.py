"""Run an analysis script inside its scientific Python environment using stdlib."""

import json
import os
import subprocess  # noqa: S404 - run the caller's explicitly supplied script
import sys
from importlib.metadata import distributions
from pathlib import Path


def main() -> None:
    """Save confirmed scientific-process facts independently of the driver.

    Raises:
        OSError: If the scientific process cannot be started.
    """
    directory = Path(sys.argv[1])
    destination = directory / "execution.json"
    command = [
        sys.executable,
        str(directory / "code.py"),
        str(directory / "inputs.json"),
        str(directory / "parameters.json"),
    ]
    record: dict[str, object] = {
        "python": sys.version,
        "interpreter": sys.executable,
        "packages": {item.metadata["Name"]: item.version for item in distributions()},
        "command": command,
    }

    def save() -> None:
        temporary = destination.with_suffix(".tmp")
        with temporary.open("w", encoding="utf-8") as output:
            json.dump(record, output, ensure_ascii=False)
            output.flush()
            os.fsync(output.fileno())
        temporary.replace(destination)

    try:
        process = subprocess.Popen(command, cwd=directory)  # noqa: S603
    except OSError as error:
        record.update(status="failed", stage="launch", diagnostic=str(error))
        save()
        raise
    record.update(status="running", pid=process.pid)
    save()
    exit_code = process.wait()
    record.update(status="completed", exit_code=exit_code)
    save()
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
