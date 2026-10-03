"""Execute native commands on an explicitly selected SSH host."""

import shlex
import subprocess  # noqa: S404 - delegate transport to the OpenSSH client
from pathlib import Path


def run_remote(
    host: str, command: tuple[str, ...], *, ssh_config: Path | None = None
) -> subprocess.CompletedProcess[bytes]:
    """Run a command over SSH without interpreting local paths or configuration.

    The remote account must use a POSIX-compatible command shell. OpenSSH owns
    authentication, host-key checking, connection configuration and exit status.

    Args:
        host: Explicit SSH destination or host alias from SSH configuration.
        command: Native argv to execute remotely, such as ``biov exec ...``.
        ssh_config: Optional native OpenSSH configuration; otherwise use its defaults.

    Returns:
        SSH result with inherited standard input, output and error.

    Raises:
        ValueError: If the destination or command is empty or malformed.
    """
    if not host or host.startswith("-") or any(c in host for c in "\x00\n\r"):
        raise ValueError("An explicit SSH destination is required")
    if not command or not command[0]:
        raise ValueError("A remote command is required")
    argv = ["ssh", "-T"]
    if ssh_config is not None:
        argv.extend(("-F", str(ssh_config.expanduser().resolve(strict=True))))
    return subprocess.run([*argv, "--", host, shlex.join(command)], check=False)  # noqa: S603
