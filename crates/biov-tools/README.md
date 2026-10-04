# Native locked tools

This crate owns a bounded local tool execution bridge, separate from CLI/MCP
adapters. It delegates dependency resolution, installation and activation to
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
receipt, failed or interrupted initial setup. A matching recorded setup is reused
without running installation again. A failed repair of an unrecorded prefix is
not a transactional upgrade and is not claimed to preserve an old working
installation. No update capability is implemented.

Setup holds an exclusive workspace operation lock; runs hold shared locks until
the manager exits. Busy operations fail clearly without waiting. These locks
coordinate this Rust API only: the Python manager and direct Pixi calls do not
participate. Roots, executable and concurrent filesystem writers must be
trusted. Paths and consistency checks are not a security sandbox. BioV does not
remove installations, caches, biological data, weights or results in this slice.
Installations are prefix-dependent, not relocatable data bundles.

CLI: `biov-rs tools setup|inspect|exec`. MCP still exposes biological datasets
and storage only; it gains no deployment-management tools. Durable lifecycle
inventory, explicit updates, uninstall and cleanup remain SPEC T56–T57 gaps.
