# Native tool adapters and locked workflows

This crate provides thin upstream lifecycle adapters and a separate bundled
locked-tool execution bridge, independent of CLI/MCP transport. It does not
implement a package solver or scientific algorithms. The initial supported tools
are `samtools` and `goatools` on Linux x86_64.

## Upstream user-tool lifecycle

CLI: `biov install`, `biov list`, `biov uninstall`. These delegate to Pixi 0.81.0
global for Samtools 1.24 and uv tool for GOATOOLS 1.6.5 with statsmodels 0.14.6.
GOATOOLS uses the complete BioV wheel's paired installed Python; a standalone
Cargo binary cannot supply that interpreter. `--uv FILE` or `BIOV_UV_BIN` selects
uv. Manager commands own environments, entrypoints, receipts and removal; no
BioV launcher, retained runner, installed registry or publication journal is added.

The dedicated layouts are `environment_root/pixi-global/{envs,bin,cache}` with an
explicit native manifest at `pixi-global/manifests/pixi-global.toml`, and
`environment_root/uv-tools/{tools,bin,cache,python}`. Always select that Pixi
manifest rather than its user/XDG fallback. These locations isolate ordinary
user-global tool settings and the separate bundled locked workspaces. The adapters
never edit profiles/global settings or change the parent shell's PATH; add the
reported dedicated bins explicitly. There is no bin-directory override.

Inventory preserves official human-readable `pixi global list` and `uv tool list`
output under root-labeled sections, consulting only existing owned roots. Manager
stdout is capped at 1 MiB and unreadable, oversized or non-UTF-8 output fails
without publishing a truncated inventory. There is no custom JSON list or
independent readiness/package-integrity audit. A first-time empty list needs no
manager and creates no root. Provisioning for locked execution does not appear in this inventory. Primary
packages are pinned; upstream managers resolve other dependencies without the
bundled scientific lock. Reinstall/retry uses their native behavior and does not
promise a fully pinned transitive set or an atomic BioV-owned transaction.

Upstream exposed commands and environments do not depend on the management BioV
wheel and survive `uv tool uninstall biov`; GOATOOLS still needs its base Python
installation to remain available. Selected uninstall removes that backend tool
prefix and its commands, preserving the other tool, caches, scientific data,
model weights, results and separate locked workspaces. Upgrades and general
cleanup are outside this slice. BioV preflights redirected backend locations and
bounded foreign/modified-command collisions; upstream entrypoints remain
authoritative and backend failures propagate without force overwrites. These
conservative guards require trusted roots and are not package authentication or
a defense against hostile coordinated native-metadata/file changes.

## Bundled locked execution

CLI: `biov tools inspect|exec`. This separate interface delegates installation
and typed JSON activation to existing Pixi 0.81.0. Both initial bundled entry
points need no preparation task.

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
trusted. Paths and consistency checks are not a security sandbox. The global
uninstall adapters never remove these locked-workflow prefixes. Software prefixes
are not relocatable biological data bundles. MCP remains a separate biological
dataset/storage surface without deployment-management tools.

Execution reads `pixi shell-hook --as-is --json`, verifies its prefix, project
root/manifest, environment name and PATH prefix, then launches the exact installed
executable with OS argv. It does not pass the executable through `pixi run`'s
command-string parsing. Empty argv stays empty even when roots contain spaces or
apostrophes. Pixi activation failures and mismatched returned paths fail closed.
Pixi evaluates activation scripts when computing its JSON environment; exported
variables reach the native executable without sourcing scripts twice.
Native `tools exec` CLI execution opts into
`execute_with_interrupt_forwarding`: direct SIGINT, SIGTERM and SIGHUP are
forwarded, and BioV reaps its native child before releasing the workspace lock.
Noninteractive execution uses a separate child process group and forwards to
ordinary descendants too. Interactive execution retains the terminal's foreground
group and forwards directly to the selected child; terminal-generated signals
retain normal group delivery. Detached descendants, signal-ignoring programs
and SIGKILL are outside this cooperative guarantee. Forwarding does not establish
exactly-once delivery across child startup or simultaneous user/group interrupts.
The general `ToolStore::execute` method preserves the embedding application's
signal policy.
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
