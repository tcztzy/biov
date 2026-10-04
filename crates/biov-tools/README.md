# Native locked tools

This crate owns a bounded local tool execution bridge, separate from CLI/MCP
adapters. It delegates dependency resolution, installation and typed JSON activation to
existing Pixi 0.81.0. It does not implement a package solver or scientific
algorithms. The initial supported bundled entry points are `samtools` and
`goatools`, on Linux x86_64 only. They do not require a preparation task.

The unchanged bundled manifest and lock are embedded in the executable and
published together into a content-keyed workspace. Its identity is SHA-256 of
manifest bytes, one NUL byte, then lock bytes, matching Python's existing
workspace identity. Cross-task reuse does not depend on the caller's analysis
directory. Existing Python workspaces can be reused after a successful native
locked setup; a native receipt is required for subsequent no-install execution.

A successful setup receipt records the exact bundled lock SHA-256, workspace
identity, manager version, platform, prefix and digest of Pixi's native prefix
marker. Inspection says `setup_recorded`, not "scientifically ready". This
records successful manager installation and checks consistency with its marker;
it does not authenticate a producer or hash/verify every installed package.
Readiness is unavailable after a changed lock, missing/changed marker, malformed
receipt, failed or interrupted initial setup. A matching recorded setup with an existing selected executable is reused
without running installation again. If that executable is absent, setup first
checks the same prefix routing under an exclusive lock, performs locked install,
and uses `pixi reinstall --locked` when normal install does not restore it.
A new receipt is published only after the selected executable is present.
Redirected prefixes and unsafe entrypoint symlinks fail instead of triggering repair. A failed repair of a damaged prefix is
not a transactional upgrade and is not claimed to preserve an old working
installation. No update capability is implemented.

Setup reuses a healthy prefix under a shared lock; installation/repair hold
an exclusive workspace operation lock; runs hold shared locks until
the selected native executable exits. Busy operations fail clearly without waiting. These locks
coordinate this Rust API only: the Python manager and direct Pixi calls do not
participate. Roots, executable and concurrent filesystem writers must be
trusted. Paths and consistency checks are not a security sandbox. BioV does not
remove installations, caches, biological data, weights or results in this slice.
Installations are prefix-dependent, not relocatable data bundles.

CLI: `biov-rs tools setup|inspect|exec`. MCP still exposes biological datasets
and storage only; it gains no deployment-management tools. Durable lifecycle
inventory, explicit updates, uninstall and cleanup remain SPEC T56–T57 gaps.

Execution reads `pixi shell-hook --as-is --json`, verifies its prefix, project
root/manifest, environment name and PATH prefix, then launches the exact installed
executable with OS argv. It does not pass the executable through `pixi run`'s
command-string parsing. Empty argv stays empty even when roots contain spaces or
apostrophes. Pixi activation failures and mismatched returned paths fail closed.
Pixi evaluates activation scripts when it computes the JSON environment; their
exported variables reach the native executable without sourcing scripts twice.
SIGINT, SIGTERM and SIGHUP sent directly to the running BioV wrapper are forwarded
and BioV reaps its native child before releasing the workspace lock. Noninteractive
execution uses a child process group, forwarding to its ordinary descendants too.
Interactive execution keeps the terminal's foreground group and forwards directly
to the selected child; terminal-generated signals retain normal group delivery.
Detached descendants, programs that ignore signals and SIGKILL are outside this
cooperative forwarding guarantee.
Forwarding is cooperative, not an exactly-once signal protocol: terminal delivery
overlapping child startup or simultaneous user/group interrupts can race.
This process-wide policy is an explicit CLI opt-in through
`execute_with_interrupt_forwarding`; the general `ToolStore::execute` method
preserves its embedding application's own signal policy.
Some punctuation-heavy roots (such as dollar-containing roots) are unsupported
by the pinned upstream activation; there is no universal pathname guarantee.

Native path arguments/environment variables are ordinary literal OS paths,
including quoted `~`; use `$HOME` or shell expansion explicitly. This differs
from the Python route's `Path.expanduser()` behavior.

This library embeds authoritative workspace assets from `src/biov/assets`.
Build from the complete checkout or sdist; standalone crates.io packaging is
not supported and this crate sets `publish = false`. The resulting installed
CLI does not need the source checkout at runtime.

Actual pinned-manager regressions are opt-in via `BIOV_TEST_REAL_PIXI`:
`cargo test --locked -p biov-tools --lib real_pixi_tests -- --ignored --nocapture`.
One uses synthetic prefixes to test activation/argv/status without downloads;
the other installs locked Samtools into a temporary root, removes its executable,
repairs it, checks unchanged lock bytes and runs the restored binary. These
complement fake-manager tests and do not establish the complete lifecycle.
