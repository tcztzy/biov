# Managed analysis: first acceptance case

This guide records the T51 contract fixed before implementation and the resulting
T52/T53 implementation and verification. The tested scope is recorded at the end.

## Scope

The first supported path runs an ordinary Python script synchronously in a
declared, locked Pixi environment on the MCP server's host. Its same-name entry
point must accept a Python script, as the `python` environment does. The existing
`biov exec` path owns installation, preparation and execution. SSH and LSF remain
available through their existing interfaces. Managed analysis does not implement
SSH dispatch or LSF submission: setting `execution_host` makes `run_analysis`
fail before creating results, rather than executing locally. This is unsupported
behavior, not just missing remote verification. There is no background scheduler
or workflow registry.

The initial real client to test is Codex CLI, with an ephemeral MCP configuration.
Official MCP SDK stdio tests are separate protocol tests and do not establish
Codex compatibility. Record the client version and actual observations when run.

The example is **complete CDS translation followed by protein properties**:

1. Explicitly select X55053.1, X62281.1, M81224.1, L31939.1 and AF297471.1 from
   Biopython's `cor6_6.gb`. Use Biopython's GenBank parser and
   `SeqFeature.translate(..., cds=True)`, genetic code 1, and check each computed
   translation against its annotation. Save all five proteins in FASTA and their
   source identifiers and parsed CDS locations in CSV.
2. Pass the complete FASTA's saved reference into a second script. Use
   `ProteinAnalysis` to compute length, average molecular weight in Da and
   theoretical pI for all five proteins. These are calculated properties, not
   experimental measurements or evidence of biological function.

The [scripts and native environment](../examples/sequence-analysis/README.md)
belong to the example workflow, not BioV's scientific algorithms.
The researcher chooses which accessions and biological questions to analyze.
Code checks exact selection, complete CDS, the declared genetic code, annotation
agreement and supported amino acids. Different accessions are different source
records, not intervals in a shared reference genome. Biopython interprets
GenBank coordinates and joined/reverse-strand features; stored parsed locations
use its 0-based, end-exclusive convention. No manual interval conversion is used.

The fixture is [Biopython 1.88 cor6_6.gb](https://raw.githubusercontent.com/biopython/biopython/biopython-188/Tests/GenBank/cor6_6.gb),
14,967 bytes, SHA-256
`01b4e193b71344752a96e2b37e886118a33c07b16c3684741c31eb3efeaf60b8`.
Its sixth record, AJ237582.1, is partial and is an explicit rejection case, not a
record silently discarded by the workflow. The example has its own minimal
native Pixi manifest and lock for macOS ARM64 and Linux x86-64; the bundled
production environments currently target Linux only.

## Request and result

One Python API, `run_analysis(AnalysisRequest(...))`, underlies the CLI
`biov analyze REQUEST.json` and MCP tool `run_analysis(request)`. The request has:

| Field | Meaning |
| --- | --- |
| `name` | Readable analysis name, 1–120 characters |
| `environment` | Declared Pixi environment, 1–80 characters |
| `code` | Complete ordinary Python script, at most 131,072 characters |
| `inputs` | At most 16 named local files or previously returned file references |
| `outputs` | 1–16 output filenames, each declared as `csv`, `fasta` or `file` |
| `parameters` | JSON parameters supplied to the script; defaults to `{}` |
| `requirements` | Inspectable input requirements and caller judgments |

Input names and output filenames are safe, flat names. They cannot name engine
files or traverse directories. Unknown request fields are rejected. This first
path accepts local files; callers can obtain identifier-backed inputs through
existing BioV access before starting it. Database access is not a prerequisite.

Each run creates a private, unique directory under `BIOV_ANALYSIS_ROOT`, defaulting
to the platform's persistent BioV data directory under `results`. It contains
`record.json`, the exact `code.py`, `inputs.json`, `parameters.json`, input copies,
execution evidence and separate stdout/stderr logs. Existing results are never
overwritten. Input files are copied, checked for changes during copying, made
read-only and hashed; scripts receive only those copies. The copies are checked
again after execution. A hash alone is not claimed to be a snapshot of a changing
source. Caller code runs with the server account's privileges; this is not a code
sandbox, and analysis scripts must not modify their inputs or prior runs.

The script receives `inputs.json` and `parameters.json` as its two arguments and
writes outputs in the run directory. An optional `checks.json` contains actual
checks as `[{"name": "...", "passed": true, "detail": "..."}]`. The example scripts
write it for both passing and failing checks. Absent scientific checks remain
unperformed; the engine never interprets exit code zero as QC success.

Responses identify `record`, `name`, `status`, `stage`, `exit_code`, `inputs`,
`outputs`, `parameters`, `requirements`, `checks`, `diagnostic`, `logs` and
`environment`. Each completed output
includes its ordinary file URI, format, byte size, SHA-256 and bounded preview.
Its execution context is this server and configured results root, not a portable
client-local path. Input references to a prior result are checked against its
completed record and digest. Missing, inaccessible or changed results fail;
neither querying nor reusing them recomputes an analysis.

The default preview shows at most two table rows or FASTA records. CSV previews
include column names and inferred preview dtypes; these are not a validation of
the whole table. FASTA previews show identifiers, lengths and sequence prefixes.
Generic files show their format/size/reference only. Order is file order, with
omissions explicit. Unknown row/record totals stay unknown. The complete response
is capped at 32 KiB, including long cells and diagnostics; full information stays
in files. This cap is independent of `BIOV_MAX_FILE_BYTES`.

`inspect_analysis(record)` reads the durable record through Python, the CLI
`biov inspect-analysis RECORD` or MCP. Querying a failed run succeeds as a query;
the original failing `run_analysis` call reports a tool execution error. MCP
protocol errors retain their native semantics. `file://` resources expose only
registered result/record/log files, with at most 1 MiB in one resource response.

For complete files larger than a response, the operator may set
`BIOV_ANALYSIS_BASE_URL` to an existing HTTP(S) storage endpoint serving the
results root. Responses then also supply `download_url`. BioV does not start a
web server, upload data or invent a chunk protocol. The client downloads from
that URL and verifies size and digest. Tests use a loopback standard-library HTTP
server serving only a temporary test directory and a separate client download
directory. No external publication is involved. A client that cannot access the
executor's files needs such a reachable endpoint; a bare executor path is not a
successful download. Retention is explicit: files remain until their owner
removes them, with no automatic cleanup or rerun overwrite.

## Execution facts and records

The driver persists launch intent before starting the existing CLI. Creating
that launcher is not proof that the scientific program started. A small stdlib
worker inside Pixi separately records the confirmed Python process and exit
status. The driver owns `record.json`; the worker owns its execution evidence.
Updates are atomic and different runs have different directories. Read-only
inspection never races a writer by rewriting the record.

The statuses distinguish launch intent, running, succeeded, failed and unknown;
failure stages distinguish launch, program and data checks. Success requires
program completion, preserved input identity, readable declared outputs and no
failed declared checks. Full original diagnostics remain in logs. Completed
outputs from the first run survive a failing second run. Partial outputs are
diagnostics, not completed data references.

After server restart, saved terminal results remain queryable. Nonterminal
records without sufficient confirmation are `unknown`; PID existence alone is
not evidence of the original process's identity or success. No query automatically
resubmits or finalizes an unconfirmed run. Cancellation or service failure can
leave execution unresolved; no execution-survival guarantee is made.

Records retain the selected environment, manifest/lock digests, actual Pixi
version, BioV driver version and interpreter, full command arguments, working
directory, exact code and parameters, inputs and outputs, and worker execution
evidence. The scientific interpreter is identified separately. BioV need not be
installed into a scientific environment merely to execute an ordinary script.
These runtime details supplement the source identifiers, genetic code, coordinate
convention and declared scientific checks.

## Acceptance cases

| Case | Required observation |
| --- | --- |
| A1: real two-step computation | Exact five selected protein IDs in input order; lengths 66, 67, 65, 65, 65; translations match annotations; second step returns five property rows despite two-row preview |
| A2: biological incompatibility | Partial CDS, missing/duplicate selections, altered translation or unsupported residues fail a named check; no silent cleaning or skipped records |
| A3: reference identity | Deleted, inaccessible or modified first-step output cannot be read or reused; no recomputation or preview substitution |
| A4: retained results | Second-step failure leaves first-step checked outputs and record unchanged; failed partial files are not completed results |
| A5: input stability | Source changes after snapshot do not affect the run; changes detected during copying or to the execution copy fail its identity check |
| A6: interruption | Interrupt before scientific launch and after confirmation; query persisted facts after server restart, never guess success/failure or resubmit |
| A7: independent runs | Repeated and concurrent invocations use different directories; neither overwrites the other's inputs, outputs or records |
| A8: bounded display | Long cells, sequence descriptions and diagnostics respect total response cap, with omissions and unknown totals explicit; changing download limits does not change preview limits |
| A9: complete retrieval | Client obtains and verifies complete bytes exceeding resource-response cap through existing HTTP storage, without using executor paths as client paths |
| A10: execution record | Required runtime fields, script, parameters, input/output identities and original logs reflect the actual run; unavailable fields are not fabricated |
| A11: MCP boundary | stdout/stderr cannot corrupt stdio; failed run has error semantics; a successful query of a failed run does not itself become an error |
| A12: real client | Codex discovers and calls analysis/query tools, sees bounded previews and can obtain full results; capture client version and distinguish this from SDK tests |

## Verification record

Verified on 2026-10-01:

| Layer | Actual observation |
| --- | --- |
| Scientific computation | Locked macOS ARM64 Pixi 0.81.0, Python 3.12.14, Biopython 1.88; five exact annotated translations and all five protein-property rows; partial CDS explicitly rejected |
| Linux execution | Docker `linux/amd64` emulation on an ARM64 host, Python 3.12.14, Pixi 0.81.0 and Biopython 1.88; all six scientific and real stdio/HTTP acceptance tests pass, including both analysis steps, restart queries, failed-step retention, invalid references and complete retrieval above 1 MiB |
| Core execution | 16 deterministic tests cover snapshots, reference identity, independent runs, checked-output retention, interruption before/after launch, post-launch record failures and bounded responses |
| MCP protocol | Real SDK stdio server processes, two-step computation, saved-record query after restart, tool-error semantics and binary resources |
| Unsupported remote configuration | An actual stdio server with `BIOV_EXECUTION_HOST` set returns a launch-stage tool error saying local execution only; no result directory is created |
| Full retrieval | SDK test downloads 1,048,832 bytes through loopback HTTP into a separate client directory; resource access above the cap fails explicitly |
| Codex client | Codex CLI 0.153.4 with ephemeral configuration calls three `run_analysis` tools and `inspect_analysis`; it reuses the full FASTA and sees two-row previews; HTTP downloads of the 272-byte properties CSV (five rows) and 1,260,000-byte synthetic transport fixture match returned sizes and SHA-256 |

The Codex run used the existing authenticated CLI without changing persistent
client configuration. Its initial sandboxed HTTP download was blocked; the same
download succeeded after network access was approved. The transport test uses
different server/client directories and HTTP, not physical machine isolation.
No claim is made that arbitrary client sandbox settings permit that access.

The macOS development driver used BioV 0.1.2 and Python 3.14.7; it is recorded
separately from the scientific interpreter. The Linux driver used BioV 0.1.2,
Python 3.12.14 and MCP 2.2.0, with all project dependencies installed from `uv.lock`.
The example's adjacent Pixi lock was used unchanged. Initial PyPI and GitHub
downloads encountered TLS connection errors. Repeating the locked dependency
installation succeeded; Pixi was installed using the existing `biov setup python
--archive` option with its official Linux archive and the shipped SHA-256 check.
Certificate verification remained enabled. The temporary container was removed
after testing; this run does not establish native ARM64 Linux or real SSH/LSF
support.

Managed SSH/LSF is not implemented. Codex desktop and other MCP clients
are not covered by these results. An attempt to inspect the Codex desktop client
was rejected by the computer-use tool's policy forbidding access to that app;
the active chat also has no BioV MCP tools loaded. Desktop acceptance therefore
remains blocked, not passed by inference from the CLI. The first path is
synchronous; it does not promise background job recovery or execution survival
after server termination.

Reproduce the independent scientific, deterministic execution and protocol tests:

```sh
uv run --locked pytest tests/test_analysis_example.py tests/test_analysis.py tests/test_analysis_mcp.py -q
```

Opt into actual Pixi/stdio/HTTP and authenticated Codex client runs separately:

```sh
BIOV_TEST_ANALYSIS=1 uv run --locked pytest tests/test_analysis_acceptance.py -q
BIOV_TEST_CODEX=1 uv run --locked pytest tests/test_analysis_codex.py -q
```

These opt-in tests require Pixi, the example environment and permission for the
temporary loopback HTTP server; the Codex test also requires authenticated model
access. The regular test suite skips them rather than presenting mock execution
as live validation.
