"""Install a child-only network filter, then exec the selected acceptance program.

Use a single-threaded launcher rather than Python preexec_fn, which is unsafe
when the pytest driver or an imported extension has already started threads.
Linux/libseccomp is required; failure to install the filter prevents execution.
"""

import ctypes
import errno
import os
import sys


def _install_network_filter() -> None:
    """Restrict only this soon-to-exec child using a Linux syscall deny filter.

    No network capability is granted and the restriction survives exec.
    Python-only monkeypatches cannot constrain the Rust binary.
    """
    library = ctypes.CDLL("libseccomp.so.2", use_errno=True)
    library.seccomp_init.argtypes = [ctypes.c_uint32]
    library.seccomp_init.restype = ctypes.c_void_p
    library.seccomp_syscall_resolve_name.argtypes = [ctypes.c_char_p]
    library.seccomp_syscall_resolve_name.restype = ctypes.c_int
    library.seccomp_rule_add.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint32,
        ctypes.c_int,
        ctypes.c_uint,
    ]
    library.seccomp_rule_add.restype = ctypes.c_int
    library.seccomp_load.argtypes = [ctypes.c_void_p]
    library.seccomp_load.restype = ctypes.c_int
    library.seccomp_release.argtypes = [ctypes.c_void_p]
    context = library.seccomp_init(0x7FFF0000)  # SCMP_ACT_ALLOW for other syscalls.
    assert context, "Cannot initialize libseccomp"
    try:
        for name in (
            "socket",
            "socketpair",
            "connect",
            "accept",
            "accept4",
            "bind",
            "listen",
            "sendto",
            "sendmsg",
            "sendmmsg",
            "recvfrom",
            "recvmsg",
            "recvmmsg",
            # Avoid native io_uring operations bypassing socket syscall rules.
            "io_uring_setup",
            "io_uring_enter",
            "io_uring_register",
        ):
            number = library.seccomp_syscall_resolve_name(name.encode("ascii"))
            if number >= 0:
                # SCMP_ACT_ERRNO(EPERM), regardless of address/family/protocol.
                assert (
                    library.seccomp_rule_add(
                        context, 0x00050000 | errno.EPERM, number, 0
                    )
                    == 0
                ), name
        assert library.seccomp_load(context) == 0, "Cannot load kernel network filter"
    finally:
        library.seccomp_release(context)


if __name__ == "__main__":
    _install_network_filter()
    os.execve(sys.argv[1], sys.argv[1:], dict(os.environ))  # noqa: S606 - Explicit test executable.
