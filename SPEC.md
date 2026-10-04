# SPEC

§I and §V describe current interfaces and contracts, including the existing
Python surface. §D describes the accepted design; §T tracks its implementation
and validation. D12's initial Rust dataset slice and paired-export reopening
extension are implemented and validated in the documented Linux source-built
scope (T63–T67). Supported scope is stated explicitly rather than inferred from
the broader design. D0 applies to every cache/data design; portable native-result
companions have their own acceptance gate (T68), with legacy migration gaps (T69).
D13 specifies the separate bounded native-source snapshot slice and its acceptance
gates (T70–T71); five-source research does not imply five implemented adapters.
D14 and T72 describe the bounded prepared RefSeq FASTA indexing extension.
D15 and T73 connect one exact prepared sequence to bounded native metric tables.

## §G GOAL
**Non-negotiable: all cache and data designs must let an agent quickly understand
and analyze the saved data without BioV.** This applies to current migration work
and every future design, including large-data storage. Standard-readable bytes
alone are insufficient if their meaning, sources or usage remain hidden in BioV.
D0 defines the invariant; current implementation gaps remain explicit.

Make BioV a biology-focused tool manager, especially for the workflows served by
`uv tool` and `uvx`: find, install, run, reuse, inspect, update and remove tools
without making users or agents remember where each task's environment was built.
Agents can already install and build software; BioV's value is owning that
software's lifecycle across tasks. Use existing package managers and scientific
software rather than replace them.

Keep biological meaning alongside this lifecycle: explicit inputs and outputs,
reference versions, coordinates, units and scientific checks make tools useful
together. Preserve complete data with bounded previews and reusable results.
Make biological identifiers a first-class entry point into reusable data, with
namespace, accession and version kept explicit. The target puts BioV's main logic
in Rust, with Rust-side Polars exposed through MCP as the primary analytical
interface. Use
Arrow IPC for table interoperability with Python and other consumers instead of
building a parallel Python analysis API. CLI and optional Python bindings expose
capabilities where useful, without requiring an agent or model. Skills
guide method selection and interpretation; repeatable handling and checks belong
in code. Task-specific orchestration and biological interpretation remain with
the agent or workflow using BioV.

## §C CONSTRAINTS
- All cache/data designs obey D0: complete standard/native data, independently discoverable meaning and lineage, portable relative references, and a tested analysis path without BioV; no exception for future large-cache designs and no claim that unmigrated caches already comply
- Target: BioV-owned main logic in Rust; Rust-side Polars/MCP analytics with Arrow IPC interoperability and optional thin Python bindings (D9–D12); no requirement for a parallel Python analysis API or legacy compatibility layers. Active development permits explicit breaking API changes. Python `>=3.12` describes the current package, not a permanent Rust-core constraint
- RuRanges-only and Biopython-only backend restrictions are superseded by contract-gated Rust migration; library choice does not change public scientific semantics
- interval semantics ≠ sequence semantics; ⊥ new DataFrame subclass or large framework
- ⊥ PyRanges objects, conversion paths, compatibility aliases, or legacy interval dependencies
- public behavior ! documented & acceptance-tested
- identifiers.org input = explicit Compact Identifier, identifiers.org URI, or BioV resource URI; bare ID inference ! curated unambiguous namespace allowlist; ⊥ registry-wide schema stripping
- identifiers.org resolution access ! fixed to `https://resolver.api.identifiers.org`; MCP transport = stdio
- Identifier-backed MCP schemes = file namespaces declared by the artifact manifest + generic `identifiers`
- Provider-data integration scales by identifier namespace × artifact kind; reuse native APIs and commands, with format conversion and input/output checks where needed; ⊥ per-function MCP wrappers; replacement of BioV-owned algorithms requires D10 validation, not a rewrite of external scientific tools
- caller-generated analysis ! use ordinary ecosystem APIs; complete script executes in selected local|LSF environment so concrete storage paths never cross executor boundary; this existing script route does not define the Rust analytical interface
- RefSeq assembly path resolution ! official `datasets download genome accession` CLI with `--include gff3,rna,cds,protein,genome,seq-report`; extracted data package structure ! preserved verbatim
- RefSeq MCP metadata ! official `datasets summary genome accession` JSON stdout; ⊥ artifact path resolution, package download/extraction or cache write
- UniProt path resolution ! official full-entry `.json` and protein `.fasta` cached independently; explicit `alphafold_cif` selects exactly one official predicted model; no implicit structure selection

## §D ACCEPTED DESIGN

D0: **Saved data must remain quickly understandable and usable for analysis by an
agent after leaving BioV.** Apply this to provider caches, derived tables,
managed-analysis outputs and any future cache/storage design. An internal
optimization may exist, but cannot be the sole usable representation or the sole
source of meaning. A live service, session handle, BioV package, repository
checkout, private database or original machine path must not be needed to read
and analyze a retained data bundle.

Keep complete data in established standard/native formats, preserving original
provider files and package layouts. Put a short readable entry point and an
inspectable machine-readable inventory beside the data, or reuse equivalent
native provider material. Identify files by paths relative to the movable bundle;
document formats, ordered schema/data dictionary, null conventions, known
scientific context, source content identities and transformation lineage. Keep
unavailable fields explicitly unknown, distinct from not applicable; do not infer
column meanings, units, coordinates, sample identity, assemblies, versions or QC
from names, dtypes or success codes. Keep declared facts distinct from verified
ones. A format/manifest schema version, BioV/tool version, biological accession
version, provider release, entry version and sequence version are separate facts.

Document an ordinary independent-reader example. Acceptance moves the bundle out
of the original location and uses standard tools with BioV absent to verify
checksums/schema, read the complete data, and filter/summarize it. Do not rely on
previews, original sources, an external/live catalog or BioV-specific wrappers.
Historical source paths can remain recorded evidence but are not portability
requirements. Source digests identify the recorded inputs; they do not imply
those inputs were bundled or independently reverified. Checksums establish
consistency against supplied metadata, not authenticity or scientific validity.

This is a design gate, not a retrospective claim of full compliance. The native
Arrow result companions in D12 are the first bounded implementation; the audited
Python/provider/fsspec/CRISPR gaps in the [migration
guide](docs/guides/rust-migration.md#data-portability-audit-and-migration-gaps)
remain open (T69). Apply the gate to each migration and every new design. It does
not expand the current 64 MiB retained-data budget or claim a large-cache engine.

D1: Users select data and an analysis, then inspect or reuse its results. The
MCP server provides analysis execution and result access using shared
core logic. Existing managed script analysis is Python-based; Rust-side
Polars/MCP is the primary analytical direction (D9, D12), with Arrow IPC for
interoperability. The client need not supply a terminal tool. Routine analysis
does not require users to choose package managers, cache paths or URI schemes.
Deployment configuration remains operator-owned; routine tool lifecycle is a
BioV product responsibility (D6–D8). Analysis execution remains distinct from
deployment management; this design does not require new MCP management tools.
V64 describes the current deployment boundary. The existing Python MCP route and
the native Rust stdio route have separately stated supported scopes.

D2: Each supported analysis declares its inputs, outputs and the metadata needed
to use them correctly. These requirements are inspectable by the caller and
distinguish fixed checks enforced in code from judgments left to the researcher.
Analysis operates on complete data and preserves complete outputs in suitable
native formats. Results have readable, unambiguous names. A result reference
identifies a saved output in its execution or storage context and can be passed
to the next analysis; preview rows are never substituted for the full input.
Missing, inaccessible or identity-mismatched results fail explicitly, without
silent recomputation, substitution or empty results. Use ordinary files and
existing storage references, without adding an analysis-specific URI scheme;
this does not exclude MCP resources using existing schemes.
Define and test how the supported client inspects, reuses and retrieves complete
results, including files exceeding a single-response limit, without assuming a
shared filesystem. Use existing storage or explicit transfer where needed; an
executor-only path is insufficient if the client cannot access it. Data may stay
on the execution host until requested; a preview need not transfer the full file.
Explicit retrieval is separate from exec's no-automatic-transfer rule (V65).
Save outputs and records at a documented location that persists across service
restarts, not solely in a disposable cache; run responses identify that location
once created. Reruns and downstream steps must not silently overwrite saved
results. Retain independently reusable outputs and records from completed steps
that passed their declared checks when a later step fails. Incomplete outputs
may be retained for diagnosis but are not presented as completed results.
Previews and records stay separate from original provider files and metadata;
existing identifier-resource contracts remain unchanged.

D3: Default analysis responses show the inputs, method, key conditions, bounded
result previews, checks performed and unresolved problems. Table previews include
column names/types, known dimensions and the first rows in the stated order;
requested filtering/sorting and omitted content are explicit. Limit both rows and
total response size, including long cells and diagnostics. Unknown totals remain
unknown rather than requiring a full scan for display. Non-tabular data retain
their native format and use appropriate summaries. Complete results and logs
remain accessible on demand; raw stdout/stderr do not become unbounded responses
or corrupt the MCP transport. Structured output alone is not a context-size limit.
Preview limits are independent of `BIOV_MAX_FILE_BYTES`, which limits downloaded
or decompressed files.

D4: Analysis readiness means satisfying the next operation's specific input
requirements, including reference version, coordinates, units and sample identity
where applicable. Coordinate conventions come from the declared native format
where unambiguous, otherwise from explicit metadata; test conversions and reject
incompatible reference versions where applicable. Skills guide method choice and
interpretation; supported format conversions and fixed checks use repeatable
code. Use reliable biological metadata and record applicable defaults; ask only
when unresolved choices affect the scientific question or interpretation.
Submission, execution completion, validation and QC outcomes are reported
separately and do not establish a scientific conclusion.
Calls may return results directly within supported client limits; work continuing
after a call returns needs a durable receipt for checking status and obtaining
results. Persist launch intent before attempting execution; record submission or
startup only after confirmation. Persist confirmed receipts and status so that
run records remain queryable after an MCP server restart; this does not imply
that execution itself survives. If interrupted before confirmation is saved,
preserve known facts and report unknown status without automatic resubmission.
Reuse protocol/executor facilities after verifying client support, without a
custom job protocol. Analysis responses expose known execution status and failure
stage (connection, launch, program or data checks) in stable structured fields,
retaining original diagnostics. Preserve MCP protocol/tool error semantics; a
successful status query may report a failed analysis. Connection loss or inability
to confirm status does not prove task failure and must not trigger automatic
resubmission of an unconfirmed run.
An inspectable analysis record retains the actual inputs and outputs, code or
command, parameters, available data/software versions and environment identity,
and logs. D0 also requires these facts to remain interpretable outside BioV.
Record input locations, byte sizes and content identities for all
sources, including local files and identifier-backed data; a path or unversioned
accession alone does not fix content. Reuse verified checksums or immutable
content versions where available; otherwise compute a digest. Establish input
stability for the supported analysis: a pre-run digest is not a snapshot or proof
that execution consumed the same bytes. Unverified consistency must be reported
and cannot count as a passed check. Output digests are required where reuse
depends on verifying unchanged bytes, rather than universally for every format.
Distinguish unavailable information from fields that do not apply; the executor
does not invent scientific metadata or infer QC success from an exit code.

D5: The first managed-script slice uses one representative two-step analysis in a configured environment,
reusing the existing Pixi manifests, locks, execution code and mature format
readers. This first locked Pixi analysis records the selected environment,
manifest and lock digests, actual execution-host Pixi and BioV versions, full
command arguments, working directory and executed code or script content where
applicable. These runtime details supplement biological metadata; Pixi-specific
fields do not become requirements for other execution methods. Inputs may be
local files or identifier-backed data; database retrieval is not a prerequisite.
The second step must consume the first step's complete output, including records
omitted from its preview. Verify real inputs and expected outputs before claiming
support. Running an arbitrary program does not establish that its results can be
interpreted or used by another analysis. Extend validated input/output handling
as needed; do not start with a workflow editor, custom query language, new
user-facing type system, persistent kernel or automatic execution-backend
switching. T51 fixes concrete MCP tool/result fields, input requirements, reference
context, each output's summary, input stability, retention and retrieval limits,
and interruption behavior in acceptance cases before T52 implementation. Reuse
existing API and format definitions where sufficient; a separate declaration
format, reference encoding or record layout is not prescribed. T51–T53 record
this implemented Python/Pixi slice; they do not require subsequent Rust table
operations to execute Python or provision a Pixi environment. D12 defines the
next bounded analytical slice.

D6: Own tool environments across tasks, not just the command that creates them.
Offer durable installation for repeated use and on-demand execution without a
persistent installation. On-demand environments may be cached and reused; they
are disposable, not necessarily deleted after every run. Expose what is installed,
where it lives, its source, resolved version/lock, platform and readiness. A tool
catalog describes available tools and their inputs/outputs; an installed inventory
describes actual local state. Reuse native manifests, locks and manager inspection
rather than introduce a second dependency solver or duplicate package catalog.
An environment is selected by its tool/dependency identity and platform, not the
analysis directory. Explicit project environments remain supported and visibly
project-owned; BioV does not silently adopt or delete them.

D7: Preserve requested exact pins and existing workflow locks; use native manager
inventory to report resolved versions and sources. User-tool installations pin
primary packages but allow upstream dependency resolution; they do not inherit
the bundled scientific workflow's full transitive lock. Reinstall/retry follows
the manager's native behavior rather than promising an unchanged resolved package
set or a BioV-owned installation transaction. Keep the bundled locked workflow's
stricter identity, receipt and execution contract separate. Reuse a compatible
installed environment before rebuilding. An unpinned request must have documented
installed/cache selection and refresh behavior,
rather than silently changing versions between tasks. Upgrades are explicit and
respect requested constraints; failed preparation or upgrade must not replace a
working installation with a partial one. Removal targets a known BioV-owned
installation. Cache cleanup targets disposable, unused content and protects active
runs and retained installations. Show locations, ownership and cleanup scope;
allow users to select storage roots and authorize deletion. Tool installations,
package/build/run caches, biological databases/model weights, and saved analysis
inputs/outputs have separate retention rules. Removing software or its cache must
not delete external data or results. Existing data-cache settings and native
manager controls remain valid; exact new controls are not yet specified.

D8: Start with a small local lifecycle slice using existing upstream managers and
one supported platform. Samtools uses isolated Pixi global install/expose/list/
uninstall; GOATOOLS uses isolated uv tool install/list/uninstall with the Python
interpreter paired with the installed BioV package. Delegate environments,
entrypoints and their native receipts/manifests to these managers. Do not create
a second installed registry, launcher/runner mechanism, publication journal or
hash-bound global version map. The independent bundled locked Pixi workflow
remains available through `tools exec`/`inspect`. Complete cross-task reuse,
inventory, explicit update, uninstall and safe cache cleanup before adding more
management layers. A tool needs one verified distribution route: community packages or official prebuilt
binaries when suitable, a container when appropriate and available, or a pinned
source build when necessary. These are alternatives, not mandatory backends for
every tool. Record the route and prerequisites; unsupported platforms or missing
build/container requirements fail clearly without silently changing sources.
Installing container runtimes, system packages or modifying shell profiles is a
separate user-controlled action. Prefer namespaced BioV entry points and native
arguments; never overwrite an unrelated executable or resolve a wrapper back to
itself. Current `exec`, `setup`, `pixi` and whole-script `run` meanings are the
starting point; D9 allows explicit interface revisions with migration guidance.
The implemented local commands are `install`, `list` and `uninstall`; future
upgrade/cleanup commands and optional short aliases remain undecided. Keep D2–D4's
biological contracts; defer broader
platform coverage, automatic backend switching and distributed environment
management until needed. The focused acceptance cases are in the
[software guide](docs/guides/environments.md#planned-tool-lifecycle).

D9: Put BioV-owned main logic in Rust. This is an explicit architecture choice
for AI-assisted engineering: use compiler-checked types, ownership and explicit
error boundaries to catch more implementation mistakes before execution. These
are reasons for the choice, not proof of scientific correctness or faster AI
development. D6–D8 borrow uv's tool-lifecycle model; Rust is not selected merely
because uv uses it. Caller scripts and external scientific tools retain their own
languages; BioV need not reimplement BWA, samtools, models or package solvers.

Treat identifiers as a distinctive, first-class capability, not strings hidden
inside fetching or analysis. An offline identifier layer owns namespace,
accession and namespace-defined version parsing, validation and canonical
representation; never split or discard accession suffixes by a global heuristic.
Do not infer biological identity or versions from a session dataset handle,
column name or arbitrary filename. Syntax validation is distinct from resolving
an identifier to a provider record or verifying that data exists. Keep existing
namespace allowlists, reserved-character handling and version-selection contracts
unless a migration explicitly revises them with acceptance evidence.

Use logical boundaries with directed dependencies:

- **identifiers**: offline namespace/accession/version types, parsing and
  validation; no provider requests, cache writes or tool execution
- **formats**: native readers/writers, schema and coordinate conversion; preserve
  native meaning, lossless representations and explicit unsupported-format errors
- **bio computation**: sequence/interval algorithms and scientific checks, using
  independently verified definitions; no transport or package-manager ownership
- **data**: provider resolution and fetching, biological datasets, provenance and
  available scientific metadata; compose identifiers, formats and cache services
- **cache**: locations, atomic storage and content/ownership/retention mechanics;
  provider selection and scientific interpretation remain outside this boundary;
  its stored data must satisfy D0 independently of BioV
- **tools**: existing manager backends, environment lifecycle and literal native
  execution, with the D6–D8 ownership and failure rules
- **optional Python**: small bindings or ecosystem adapters when justified;
  no independent analytical implementation or mandatory compatibility facade
- **CLI/MCP**: thin command/protocol adapters over Rust operations; no duplicated
  scientific, identifier, provider or lifecycle rules

These are responsibility boundaries, not an instruction to create one crate per
bullet. Use modules or crates according to actual dependencies and substantial
responsibilities; there is no fixed tiny crate budget and no requirement to
create empty crates in advance. The initial core/PyO3 workspace is a starting
point, not a permanent architectural limit.

Rust-side Polars is the chosen primary engine for columnar analysis through MCP.
Operations consume complete datasets, return reusable dataset references and
provide schema, row counts when known, and bounded previews. Arrow IPC is the
interoperability boundary for complete typed tables, including consumers using
Python Polars, PyArrow or another Arrow reader. Do not build a second full Python
analysis API or require Python to call the Rust MCP implementation. The existing
tiny PyO3 sequence binding may remain useful without expanding into a parallel
analysis surface. Do not emulate pandas indexes, extension dtypes or subclass
behavior solely for compatibility.

Select schemas, sequence representation, null/NaN handling, supported operators
and row-order rules explicitly for each slice. Conversion must preserve values,
duplicates and scientifically relevant metadata, or reject unsupported inputs;
never silently stringify or drop them. Arrow IPC alone does not establish
reference versions, coordinates, units or sample identity. Retain known metadata
in a documented record or representation and mark unavailable metadata as unknown.
Rust-Bio, noodles and Rust interval kernels remain candidates, not assumed
drop-in replacements. Check behavior, license, maintenance, platform/build cost
and end-to-end performance; neither Arrow nor Rust implies zero-copy or a speedup.

Existing Python usage informs scientific requirements without freezing names,
signatures or return types. Explicit breaking updates are acceptable. Keep the
legacy Python `biov mcp` route separate while the new `biov mcp-native --data-root
DIR --output-root DIR` route is implemented and validated with the official Rust
MCP SDK (`rmcp`) over stdio. Do not silently reroute the existing entry point or
claim full Python-server feature parity, Python-free BioV distribution, or
client acceptance from a successful compile. No new per-operation RPC framework
or universal backend abstraction is required. See the [migration
guide](docs/guides/rust-migration.md) for actual interfaces and validation scope.

D10: Migrate by scientific and operational contract, not by a blanket rule to
keep old dependencies or a single wholesale rewrite. Inventory relevant current
outputs, schemas, metadata, warnings/errors and CLI/MCP behavior, including Python
types only where that slice still exposes them. Specify the retained or
intentionally revised contract in offline fixtures; a new Rust/MCP slice need
not first reproduce every legacy Python API. Pin the reference
implementation/version and record reference assemblies,
coordinates, units, genetic codes and numerical tolerances. Preserve the
scientific meaning in V1–V13 (including coordinates, weighted GC and protein
units), identifier/version handling, original provider bytes, executor
locality, interruption semantics and D2–D4 result integrity. Compare complete
results, never just bounded previews. Old behavior is a differential reference,
not an API freeze or proof that a scientific result is correct: independently
check published definitions, authoritative tables and hand-calculated or
independently implemented examples. Resolve discovered bugs as explicit
contract revisions rather than silently reproducing or fixing them.

Replace sequence algorithms with existing Rust libraries where verified. A
small BioV-owned implementation is allowed for gaps, but needs documented
mathematics, source/provenance for constants and tables, independent oracle
cases, adversarial and property tests and review before it replaces the
reference. AI-generated code has exactly the same gate. Rust memory safety and
compilation do not validate biology. Biopython-backed `Seq`/`SeqRecord` return
types may be replaced outright with a documented Rust/MCP representation and
Arrow IPC interchange where appropriate; no Python mirror or compatibility shim
is required. Remove the dependency when no remaining BioV function needs it. Scientific Python dependencies in caller-owned analysis
environments need not disappear when BioV's core moves to Rust.

D11: Validate and distribute the Rust binary independently of optional Python
bindings. Pin the Rust toolchain/MSRV, selected dependency features and Cargo
lockfile; test Rust operations without Python. The native CLI/MCP release gate
includes argv/help/exit/output contracts, clean stdio, protocol initialization,
tool discovery, structured results and errors, and complete-data reuse/export
through a real MCP client. The official `rmcp` SDK supplies protocol machinery;
a custom JSON-RPC implementation is not the product. A tested narrow native
analytical server does not establish parity with the existing provider,
identifier or managed-Python-analysis server.

The complete Python distribution uses upstream setuptools-rust to build both
`biov._native` and the single native `biov` executable. It has no Python `biov`
console-script owner and no separate `biov-rs` executable. Maturin's PyO3 target
selection does not package the binary alongside this extension; do not hand-edit
wheel RECORD files or inject binaries with a private packaging backend.
Existing Python capabilities are explicitly delegated to the installed package's
paired interpreter. Native routes never start Python.
Where bindings are shipped, test installed-wheel imports and selected Python/CLI
contracts; Python 3.12–3.14 is the current acceptance matrix, not a dependency of
the standalone Rust analytical server. Declare and test any revised Python
minimum before release. Evaluate stable ABI against actual binding/dependency
needs; do not assume abi3 covers every runtime, free-threaded build or OS/CPU.
Do not expand the binding solely to keep an old Python analysis API intact.

Initial native release targets are Linux x86_64 and macOS ARM64. Each binary or
wheel requires actual install/runtime tests for its target, with a stated Linux
libc baseline and macOS deployment target. Other platforms remain unclaimed
until validated. Existing users outside that set need a documented source-build
route before the native release. Supported prebuilt-artifact users need no Rust
compiler; source/sdist builders need the stated Rust/linker/system prerequisites.
For Python installations, check retained dependencies' wheel availability as well
as BioV's before claiming compiler-free installation.

Include all required Rust sources, Cargo metadata/lock, Python wrappers,
assets/licenses, tests and docs in a self-contained Python sdist, and build/test
from the unpacked sdist. Standalone CLI archives are separate artifacts, with
checksums and their own target/dependency tests; wheels are not universally
portable binaries. Keep scientific environment platform support separate from
BioV package support. Do not advertise Rust speedups without reproducible
release-build, representative-data measurements of conversion, memory and runtime.

D12: The implemented analytical slice is a local Rust/Polars dataset workflow
over `biov mcp-native --data-root DIR --output-root DIR`, using the official `rmcp`
stdio server. The initial CSV/query/export scope is validated on Linux x86_64,
including an installed binary outside the checkout (T63–T65). Paired-export
reopening is implemented and source-built Linux acceptance covers two real MCP
processes (T66–T67). Desktop-client integration and portable prebuilt releases
remain unclaimed. The seven tools are `dataset_open`, `dataset_reopen`,
`dataset_preview`, `dataset_query`, `dataset_export`, `dataset_read_artifact`
and `dataset_release`. The server opens local UTF-8 CSV files with a header,
provides schema, exact materialized row count and bounded previews, creates
derived datasets from complete-data operations, exports Arrow IPC and reopens
previously exported JSON records with their paired IPC files. No arbitrary Arrow
IPC import, TSV, Parquet, aggregation, arbitrary SQL or code execution is claimed
in this slice. This is a data-handling foundation, not biological
normalization, sequence/interval parity or completed tool lifecycle management.

Opening reads one bounded in-memory byte snapshot and computes the input SHA-256
from those same bytes before parsing. `dataset_open` accepts an optional `schema`
map from existing column names to `string`, `int64`, `float64` or `boolean`.
Unspecified columns remain strings: never infer numeric types and lose leading
zeros in identifiers or the exact spelling/precision of numeric-looking text.
Numeric and Boolean operations require an explicit column declaration. Invalid
values for the declared type and non-finite Float64 values fail explicitly;
preserve such source text as strings or clean it in an explicit separate step.
Preserve nulls rather than converting them into empty strings or NaN.

Before allocating Polars columns, the CSV preflight validates UTF-8, header names,
column count, duplicate headers, declared-schema names, every row's field count,
declared values and a conservative payload/cell allocation bound against the
remaining session budget. One established `csv-core` decoder determines fields
and records in both preflight and direct typed construction; Polars does not
reparse CSV bytes. LF, CRLF, CR and mixed record terminators are accepted. Empty
physical lines outside quotes are ignored; quoted embedded line endings and
blank lines retain their exact bytes. Doubled quotes decode to one quote;
unterminated quoted fields and characters after their closing quote fail.
Unquoted empty fields are null; quoted empty strings stay distinct from null.
Empty numeric/boolean fields are null even when quoted. Reject ragged rows
rather than padding or truncating them, and verify constructed cardinality
against the preflight plan. The complete dialect is documented in
`docs/guides/rust-datasets.md`. One query applies an optional typed filter (`eq`,
`gt`, `ge`, `lt`, `le` or `is_null`), then one optional stable sort with nulls
last, then optional column selection. Boolean filters support only `eq` and `is_null`; other comparisons
require values of the selected column's supported type. Unknown/duplicate
selected columns and unsupported operators fail explicitly. No implicit row
limit is applied to query results. Declare any revision of these rules before
claiming it as supported.

Dataset handles identify immutable logical inputs/results within one server
session. They are neither biological identifiers nor durable saved-result
references; unknown, released or expired handles fail explicitly. Original and
derived handles remain separate, and filtering never operates on preview rows.
The first offline Rust identifier subset validates explicit RefSeq GCF and
UniProt references, with bare GCF as the only inference route. It is not a port
of the whole identifiers.org registry or prompt scanner and does not resolve
providers. Namespace-specific versions remain explicit; a UniProt accession
alone does not establish its entry or sequence version.

Caller-declared metadata fields are `identifier`, `species`, `reference`,
`coordinates` and `units`. Missing values are null/unknown, and supplied values
are declared rather than provider-verified. A supplied identifier receives local
syntax validation, not evidence of record existence or content identity. Do not
infer scientific context or sample identity from CSV content, propagate metadata
as if provider-verified, or report unperformed biological checks as passed.

The concrete initial resource limits are:

- CSV input: 16 MiB; one to 64 columns, each name one to 128 UTF-8 bytes
- Paired reopen: 64 KiB JSON record and 65 MiB IPC byte snapshot; export checks
  these same byte caps before publishing either file; at most 1 MiB
  per IPC footer/message metadata block, 4096 record batches and 4096 buffers
  per batch. Record/file paths are at most 1024 bytes. These are parser limits,
  not an increase in the retained-data budget or a process-memory ceiling
- Session datasets: at most 16 handles and a 64 MiB retained-data charge, summing
  each Polars estimated size plus 17 bytes per cell to cover omitted string-view
  buffers and validity overhead conservatively. Before Polars allocation, CSV
  preflight combines decoded string bytes, numeric payload, packed boolean and
  validity bounds, this per-cell charge and already-retained datasets. Temporary
  source/parser buffers are outside the retained charge. IPC preflight combines decoded
  logical payload, 17 bytes per cell and the currently retained charge before
  full decoding. This is not a peak-memory ceiling
- Preview: zero to 50 rows, default five; 24 KiB serialized-row budget and
  256-byte text-cell limit, with omitted rows and truncated cells explicit
- Derivation: at most 32 recorded steps and 8 KiB serialized provenance on open and query;
  metadata values at most 256 bytes each
- Exports: at most 64 registered artifacts per session; UUID-generated filenames
  and create-new publication, never a caller-controlled arbitrary destination
- Explicit artifact reads: one to 48 KiB raw bytes per call, base64-encoded with
  offsets, total size, SHA-256 and end-of-file status

Preview rows alone do not bound a whole response. The adapter additionally caps
serialized tool results at 95 KiB, text at 16 KiB and errors at 2 KiB; tests
verify normal JSON-RPC envelopes below 96 KiB, including maximum-size artifact
chunks. Recorded parameter values are constrained by the provenance cap. A failed operation must not
be registered or reported as a completed result. Complete files left after a
later record-publication failure may remain for diagnosis, without a success
claim or implicit retry.

Roots must already exist, have UTF-8 paths and be trusted, operator-owned
directories. Input paths are data-root-relative; reject absolute paths, traversal
and symlink escapes. Generated exports remain under the configured output root
and do not mutate
source files or overwrite earlier exports. Canonical path checks establish the
supported boundary, not a secure sandbox against hostile concurrent local
filesystem writers. Do not offer arbitrary input URLs or caller-supplied export
paths through this surface.

Each export saves the complete Arrow IPC file and an adjacent JSON record with
schema, row count, source snapshot identity, transformation provenance, declared
metadata, software information, output byte count and SHA-256. The upstream
Polars IPC writer's oldest compatibility setting emits strings as Arrow
`LargeUtf8` (`large_string` in PyArrow), supporting direct ordinary-reader
filtering without a BioV wrapper or reader-side string cast. Other primitive
logical types and strict CSV record v2 fields are unchanged. Files persist
independently of the server session and are not automatically deleted on dataset
release or shutdown. Dataset and artifact handles both expire with the session.
`dataset_release` frees only a dataset handle and leaves saved exports intact.

New exports additionally save `artifact_<id>.manifest.json` (manifest version 1)
and `artifact_<id>.README.md` beside the existing Arrow and strict version-2
CSV record. The manifest uses same-directory basenames for Arrow/record/README,
records Arrow and record byte sizes/SHA-256, complete row count, ordered
column types/null counts and explicitly unknown column semantics. It retains
caller-declared scientific metadata, syntax-only parsed identifier/version facts,
source identity and historical location, ordered operations, declared CSV schema,
software facts and reopening context. Source paths are historical references,
not required bundle members; the original source is not implicitly included.
Unknown meanings, column units/coordinates and independent reference/provider/
entry/sequence versions remain null. Format/manifest versions and software
versions never fill those biological fields. See the [native dataset
guide](docs/guides/rust-datasets.md#use-a-result-without-biov) for field details.

The README has a standard-PyArrow example requiring no BioV: verify file sizes
and hashes, read all rows, check schema/count, filter and summarize. These
companions are descriptive, untrusted data; paired reopening still checks only
the strict JSON record and IPC and does not require or validate companions.
The existing strict CSV record v2 schema is unchanged; D15 separately introduces
record v3 for sequence metrics. Companion limits are 128 KiB
for the manifest and 16 KiB for the README, checked with the existing caps before
any final filename is published. Publication of all four files is not a single
atomic transaction; I/O failure can leave a partial group, but no completed
export or artifact handle is reported before all four are saved. The MCP export
response adds host `manifest_path`/`readme_path`, not their complete bodies;
`dataset_read_artifact` remains Arrow-only. A remote client needs an authorized
ordinary file transfer to acquire the companions. Host paths alone do not make
a complete bundle available remotely.

The paired-export reopening contract is `dataset_reopen({record_path, preview_rows})`,
where `record_path` is a `.json` record path relative to the configured data root and
`preview_rows` defaults to five with the normal preview bounds. Reopening requires
both files: the record's `file` must be a basename naming the IPC file in the same
directory, matching the native `artifact_<32 lowercase hexadecimal digits>.arrow`
name and recorded artifact ID. An in-root record symlink uses the canonical
record's parent directory to locate its IPC file. Both resolved paths must remain
inside the data root. Reject
absolute paths, parent traversal, symlink escapes and arbitrary IPC-only input.
Read the bounded IPC bytes once; compute the SHA-256 and parse the complete table
from that same in-memory snapshot. Validate record version and format, byte size,
SHA-256, ordered schema against actual supported column types, row count and
metadata fields before publishing a new session handle. Missing, malformed,
unsupported or mismatched pairs fail explicitly, without reparsing a CSV source,
silently recomputing the result or importing unchecked IPC. Retained dataset,
preview and provenance budgets still apply. IPC is restricted to native
little-endian, uncompressed, flat string/int64/float64/boolean exports; preflight
rejects dictionaries, extensions, other types and unsupported metadata before
the full Polars schema/array allocation. String storage accepts both current
`LargeUtf8` and the prior native `Utf8View` encoding. This remains a narrow
preflight using upstream Arrow metadata definitions before the standard Polars
reader, not a new generic Arrow import implementation. For `LargeUtf8`, require
exactly `(rows + 1) * 8` offset bytes, a zero first offset, nonnegative monotonic
offsets within the values buffer, valid UTF-8 with every offset on a character
boundary, and a final offset equal to the values-buffer length. Charge complete
string payload plus the existing 17 bytes/cell and retained datasets under the
same 64 MiB budget; no bound is relaxed for compatibility.

Read initial-slice CSV record version 1 and version 2; export CSV version 2 with a
`reopen_verification` field that is null for CSV-opened datasets. Re-exporting
old native `Utf8View` uses the upstream compatible writer and emits `LargeUtf8`.
This is one-way runtime compatibility: earlier development readers supporting
only `Utf8View` cannot reopen these newer `LargeUtf8`
exports. The unchanged record schema version describes logical record fields,
not a guarantee that an older binary accepts a newly supported physical encoding.

Successful reopening exposes `reopen_verification` and carries its verification
context through queries into the next export record. The context includes the
record path and SHA-256, artifact SHA-256 and bytes, record version, the checks
`artifact_bytes`, `artifact_sha256`, `schema` and `row_count`, and the states
`input_consistency: parsed_same_in_memory_snapshot_as_digest`,
`original_provenance: recorded_claims_not_independently_verified` and
`authenticity: not_established`. A later reopen replaces this context with its
current checks instead of nesting prior verification objects. These checks establish
consistency between the supplied record and the exact table bytes consumed.
The record is not a trusted signature or security attestation: changing the IPC
and updating its record together can pass these checks. Recorded source history,
software claims and biological metadata remain recorded/unverified; file
integrity does not establish provider authenticity, biological correctness or
QC. Do not silently promote caller declarations into verified scientific facts.

After a restart, configure the old export directory as the new data root and
explicitly reopen a saved pair to obtain a fresh handle. Reopening neither
resurrects old dataset/artifact handles nor reconstructs a persistent catalog.
The saved JSON and IPC files are the reusable artifacts; no daemon database,
background job, automatic recovery or general-purpose import mode is introduced.

An execution-host path is not proof that a remote MCP client has the file.
`dataset_read_artifact` provides explicit bounded transfer during the originating
session, verifying the saved bytes against their export identity before returning
a chunk. Validate full client-side reconstruction and independent Arrow readback,
including rows omitted from previews. Reopening a pair already present on the
execution host does not transfer it from a remote client or restore an expired
artifact retrieval handle; export the reopened or queried dataset to obtain a
new session artifact handle. Broader provider ingestion, durable catalogues,
automatic restart-time retrieval, scientific metadata inference, lifecycle
additions and retirement of the Python route require separate contracts and validation.

### D12 acceptance cases

The initial local-table cases are validated by T63–T65. T66–T67 validate
paired-export reopening and two-process reuse in the Linux source-built scope.
T68 separately gates the portable companion bundle with no BioV in the reader.

```gherkin
Feature: Complete local datasets through Rust MCP and Arrow IPC
  Background:
    Given biov is running as an rmcp stdio server
    And a data root and a separate output root are configured

  Scenario: Filter complete data and independently read back its Arrow export
    Given a CSV with columns record_id and score contains these ordered rows
      | record_id | score |
      | first     | 1     |
      | second    | 2     |
      | retained  | 9     |
      | retained  | 9     |
      | last      | 7     |
    When I open the CSV with score declared int64 and request two preview rows
    Then the schema and total row count of five are reported
    And the preview contains only first and second with omissions declared
    When I filter the dataset for score greater than 5
    And sort score descending then select record_id and score
    Then a separate derived dataset has three rows
    And the original dataset still has five rows
    When I export the derived dataset to Arrow IPC under the output root
    And retrieve every bounded artifact chunk through MCP
    And reconstruct the complete file and verify its declared SHA-256
    And read that file with an independent Arrow-compatible reader
    Then all three rows appear in the requested order
    And both duplicate retained rows are preserved
    And the complete score column is 9, 9 and 7
    And the selected column names, types and values match the complete result
    And a persistent JSON record identifies the source snapshot and operations
    When I release the dataset handle
    Then the handle is unusable and the exported files remain intact

  Scenario: Analyze a moved result bundle without BioV
    Given an exported Arrow, strict record, descriptive manifest and README
    And the table contains typed values, nulls, leading-zero IDs and duplicate rows
    When I move those four files to a new directory and remove the original source
    And I use the README's ordinary-reader example with BioV absent
    Then only bundle-relative filenames are needed to find the result files
    And Arrow and strict-record byte counts and checksums agree with the manifest
    And the complete table schema, rows, types, nulls, duplicates and order match
    And filtering and numeric summaries include rows absent from the original preview
    And known declarations and source hashes are inspectable without BioV
    And unavailable meanings, scientific context and biological versions remain unknown
    And consistency checks do not authenticate the producer or validate biology

  Scenario: Reopen a complete saved result in a second MCP process
    Given MCP process A opens a CSV with explicitly declared primitive types
    And the complete input contains identifiers with leading zeros, duplicate rows and nulls
    When A filters, stably sorts and selects columns from the complete dataset
    And A exports the derived result as an IPC file and adjacent JSON record
    And I reconstruct every exported IPC byte and independently read the complete table
    And I close A's stdin and wait for its clean exit
    And I start MCP process B with A's output directory as its data root and a new output directory
    Then A's dataset and artifact handles cannot be used in B
    When B calls dataset_reopen with the saved record's relative path and preview_rows 1
    Then B returns a fresh dataset handle, exact schema and full row count
    And reopen_verification describes record-to-byte and record-to-table consistency checks
    And the preview does not replace any rows in the reopened dataset
    When B previews, queries and exports the reopened dataset
    And I reconstruct B's full export, verify its SHA-256 and independently read it
    Then every expected row, column type, null, duplicate and ordering matches
    And the new export record preserves the reopening verification context
    And recorded provenance and biological metadata are not reported as authenticated

  Scenario Outline: Refuse an unusable export pair without fallback
    Given a saved JSON record and its adjacent IPC file
    And the pair <problem>
    When I call dataset_reopen with the record's data-root-relative path
    Then I receive a bounded explicit tool error and no new dataset handle
    And no original input is reloaded, result recomputed or unchecked IPC imported
    And existing completed exports remain unchanged

    Examples:
      | problem                                                 |
      | has a missing record or IPC file                         |
      | contains malformed JSON or an unsupported record version |
      | declares an unsupported format                          |
      | names an absolute or non-basename IPC path               |
      | resolves either file outside the data root               |
      | exceeds a documented record or IPC limit                 |
      | has a byte size or SHA-256 mismatch                       |
      | has a column name, order or actual type mismatch          |
      | has a row count mismatch                                 |
      | has invalid scientific metadata fields                   |

  Scenario: Preserve identifier spelling and numeric-looking text by default
    Given a CSV row has record_id 001 and measurement 9007199254740993.0010
    When I open it without declaring column types
    Then both columns have string type
    And both values retain their complete original spelling
    When I export and independently read back the complete Arrow IPC file
    Then the values are still 001 and 9007199254740993.0010 as strings
    And numeric comparison requires an explicitly declared numeric column

  Scenario: Do not invent biological identity or metadata
    Given a local CSV with no declared biological metadata
    When I open it and derive a filtered dataset
    Then its dataset handles are session references rather than biological IDs
    And reference assembly, coordinates, units and sample identity remain unknown
    And no biological validation is reported as passed
    When I use an unknown or expired dataset handle
    Then I receive a clear tool error and no substitute dataset

  Scenario Outline: Reject limits and unsafe paths without publishing success
    Given a request that <violation>
    When I call the relevant dataset operation
    Then I receive a bounded explicit tool error
    And no result is registered or reported as completed by the failed operation
    And existing inputs and completed results remain unchanged

    Examples:
      | violation                                      |
      | exceeds the documented input-byte limit        |
      | exceeds the documented retained-dataset limit  |
      | requests a preview beyond the documented limit |
      | has duplicate CSV headers                      |
      | has a ragged CSV row                            |
      | declares a schema for an unknown CSV column     |
      | has an invalid value for its declared type      |
      | contains a non-finite declared float64 value     |
      | exceeds the conservative cell-allocation limit  |
      | exceeds the 8 KiB provenance limit         |
      | selects an unknown column                      |
      | uses an unsupported filter operation           |
      | opens a path outside the data root             |
      | opens a symlink escaping the data root         |
      | supplies an arbitrary export path              |
      | exceeds the registered-artifact limit          |
      | requests more than 48 KiB in one artifact read  |
      | reads an artifact modified since export        |
```

D13: Native storage preserves complete provider packages as immutable, independently
usable filesystem snapshots. The initial `biov-storage` Rust library copies an
already-present package from an explicit trusted source root into a separate
trusted store root. It never moves, deletes, hardlinks or rewrites original source
files. Native metadata/layouts remain authoritative; a small operational receipt,
relative checksum inventory and wrapper README add only missing navigation,
identity and registration facts. The README directly lists analysis entry files,
representations and scope so the next reader need not reverse engineer the package.
The complete data must remain analyzable without BioV, its database, an original
source path or network. See the [native-storage guide](docs/guides/native-storage.md).

The first semantic adapter is materialized NCBI Datasets RefSeq: validate an exact
canonical versioned GCF against the native catalog, allow accessionless report
groups, resolve catalog paths relative to `ncbi_dataset/data`, require native
README/checksums/catalog, verify catalog member sizes and provider MD5 entries,
and expose only actual available representations. This does not validate all
biological relationships in FASTA/GFF. PDB registration is explicitly
caller-declared identity, entry or assembly scope and representation/path mappings;
syntax, paths and bytes are checked, while mmCIF semantics and claimed biological
identity are not. UniProt, AlphaFold and GEO semantic adapters remain planned.

Native-source inventory identity includes the complete sorted relative file and
directory tree, sizes and file SHA-256 values using a documented canonical
encoding. Wrapper files do not affect that identity. The readable layout is
`artifacts/<namespace>/<canonical accession>/snapshots/sha256-<full digest>/` with
`source/`, `acquisition.json`, `checksums.sha256` and `README.md`. Receipt version,
registration time, biological accession version and provider version/release facts
are distinct. Unknown original acquisition URL, time and client stay null; no
historical source location is required. Conflicting declarations for identical
source content fail without mutating the saved receipt.

Copy and verify in staging, then publish atomically without replacing a winner.
Use the documented local locking and no-replace filesystem primitive. Concurrent
identical copies reuse a verified existing snapshot; changed annotation under the
same assembly accession creates a new snapshot. Interrupted/partial stages are
never ready. An unsuccessful registration leaves the prior valid snapshot intact.
This is a trusted local-filesystem boundary, not protection from hostile writers
or proof of multi-host/distributed publication semantics.

Offline resolution scans the filesystem, checks recorded source bytes and returns
ordinary host paths plus portable store-relative paths and bounded snapshot
summaries. No DNS or download is attempted. Missing identity, multiple matching
snapshots, unavailable representation and corrupted data are explicit outcomes.
Selection is never implicit newest-by-mtime; callers pin an exact snapshot when
more than one matches. Moving the store must allow fresh discovery with no
original source root. No durable discovery index exists in this slice, so this
is reconstructed filesystem discovery, not a database rebuild feature.

Source copies/hashes use bounded chunks, but metadata/traversal/output cardinality
have explicit limits. Full verification on resolve may be expensive. This does
not expand D12's independent 64 MiB retained-data charge, guarantee peak memory,
provide lazy analytical scans, or establish large-table/distributed support.
Downloads, archive acquisition/hydration, automatic cache/import migration,
external registration, aliases, GC, quotas and a persistent index remain separate
future work. Original Python artifact paths retain their existing behavior.

An analysis-ready derived layer is planned separately from preserved native
sources. It should expose useful complete tables, indices and relationships with
input hashes, actual commands/software, schemas, units/coordinates and unknowns.
The planned layers are acquired native data, prepared/materialized analysis-ready
views and durable results. Transform identity includes input hashes, code/tool
versions, parameters, seeds and output-affecting dependencies; output SHA-256 is
separate. Nondeterminism or uncaptured external state prevents deterministic-reuse
claims. Prepared views need row meaning, joins/shard order and null semantics.
Examples include RefSeq genome/annotation indices, GEO expression/sample/probe
views retaining MAS5 meaning and multi-mapping, and PDB entry/assembly views.
Reuse immutable raw references without mandatory payload duplication. Portable
export must explicitly materialize dependencies or state its output-only scope;
current registration does not convert formats or rewrite custom metadata. Future
metadata/scan/indexed/materialize interfaces require a memory budget and moved-
bundle/cache-invalidation acceptance; Arrow is not a universal conversion target
for FASTA/BAM or other domain-native representations.

Preservation is separate from decoding: a future generic file catalog may keep
opaque formats with explicit reader capabilities or unsupported-decoding status.
Never automatically unpickle or execute untrusted stored content. A successful
registration or `ready` file resolution proves the declared local materialization
and recorded integrity checks, not format parseability or scientific validity.

### D13 acceptance cases

These are the contract gates for the bounded slice. The native-storage guide and
T70 record actual executed coverage; real provider inspections are evidence for
the examples, not proof of every lifecycle or scientific edge case.

```gherkin
Feature: Offline source-native snapshots usable without BioV
  Background:
    Given an existing trusted source root and a disjoint trusted store root
    And source originals remain outside the managed snapshot tree

  Scenario: Register and resolve a complete native RefSeq package
    Given a materialized NCBI Datasets package for GCF_000005845.2
    And its native README, MD5 inventory, catalog and catalog members agree
    When I register the package with its exact canonical versioned reference
    Then the complete source tree and bytes are copied without rewriting
    And the wrapper directly lists analysis entry files and their representations
    And original acquisition facts not known from registration remain unknown
    When I resolve genome_fasta without network access
    Then I receive an ordinary local path to its native catalog-selected file
    And the source original is unchanged

  Scenario: Distinguish representation availability from biological absence
    Given a valid RefSeq snapshot with no RNA FASTA entry in its native catalog
    When I resolve rna_fasta offline
    Then the result is unavailable with the representations actually present
    And no download occurs and no biological absence is inferred

  Scenario: Concurrent identical registrations cannot replace a winner
    Given two writers registering the same complete native source tree
    When both attempt to publish the same snapshot identity
    Then one complete verified snapshot is retained
    And the other writer reuses the verified winner
    And its immutable receipt and source bytes are not replaced

  Scenario: Changed annotation does not overwrite the same assembly version
    Given a verified snapshot for GCF_000005845.2
    When I register a package with changed annotation bytes for GCF_000005845.2
    And its native checksums and catalog reflect the changed package
    Then a distinct content-addressed snapshot is published
    And the old snapshot remains byte-identical
    When I resolve without a snapshot selector
    Then I receive an explicit ambiguity with bounded choices
    When I select the original snapshot exactly
    Then its original native representation remains available

  Scenario: Partial staging and invalid input do not become ready
    Given a valid snapshot and an interrupted partial directory in staging
    When I rediscover the store and attempt to register an invalid package
    Then staging is ignored as a resolution candidate
    And the invalid registration does not publish a ready snapshot
    And the earlier valid snapshot remains intact

  Scenario: Reconstruct discovery after moving the store
    Given a published store with no durable discovery database
    When I move the complete store to another directory
    And the original source and store locations are no longer available
    And I open a new store instance at the new root
    Then scanning reconstructs the same snapshot identities and representations
    And ready paths point inside the moved store
    And no original path, old session or network lookup is needed

  Scenario: Read and analyze a copied snapshot with an ordinary reader
    Given a complete copied snapshot outside its original store
    When I use standard readers in an environment without BioV and with network blocked
    Then relative checksum paths verify the complete native file bytes
    And the native metadata explains the representations and scientific scope
    And I read complete records and perform a meaningful summary or filter
    And no original source path or private catalog is required

  Scenario: Keep PDB entry and assembly claims explicit
    Given a caller-declared PDB package scoped to entry coordinates
    When I register and resolve its explicitly mapped coordinate representation
    Then the ordinary native file is returned with entry scope
    And validation states that biological identity and mmCIF semantics are unverified
    And it is not silently selected for an assembly-specific request

  Scenario: Reject corruption and declaration conflicts without repair by overwrite
    Given a previously published snapshot
    When retained source bytes no longer match its receipt
    Then resolution reports corrupt rather than a ready path
    And a new registration does not silently overwrite that snapshot
    When identical native source bytes are registered with conflicting declarations
    Then a declaration conflict is reported without changing the existing meaning
```

D14: Prepared RefSeq FASTA indexing is a separate Rust library responsibility,
consuming D13 exact verified snapshots. It accepts a canonical versioned RefSeq
reference, exact source snapshot ID and an explicitly selected `genome_fasta`
source path. Plain uncompressed case-preserving IUPAC DNA FASTA is indexed with
pinned upstream noodles-fasta; BioV adds bounded input validation, duplicate-ID
rejection, provenance and publication, not a replacement FASTA/index decoder.
The immutable native package is unchanged. A conventional FAI index, small TSV
sequence dictionary, readable README and JSON provenance are published below
`prepared/sha256-<recipe>/`. Schema, identifier token semantics, base versus byte
units, interval conventions, unsupported input and unknown scientific facts are
explicit. No full sequence table or mandatory raw-data duplication is introduced.

The recipe identity includes exact selected input bytes/snapshot/path,
implementation contract and output-affecting dependency versions/parameters.
Actual output SHA-256 values remain separate. Output-affecting algorithm changes
must update the implementation contract. Publication is staged, verified and
atomic without replacement. Reuse verifies sources and outputs, reports corruption
rather than overwriting, and never treats incomplete staging as a ready artifact.
Sequential validation has explicit line, metadata and record bounds; it does not
claim constant-time reuse, whole-process peak memory or benchmarked performance.

Relative references reuse the native snapshot. A portable dependency closure
contains the exact referenced complete snapshot and prepared files at their
store-relative paths. The prepared directory alone is explicitly incomplete.
Ordinary indexed readers must work after that closure moves, without BioV, a
catalog database, original source locations or network. Acquisition facts remain
native/recorded facts; hashes establish consistency, not authenticity or QC.
General transforms, compressed FASTA, GEO conversion, downloads and GC remain
separately planned. See `docs/guides/prepared-fasta.md` for the precise implemented
contract and executed validation status.

### D14 acceptance cases

```gherkin
Feature: Independently readable prepared reference indexing
  Scenario: Prepare, reuse and read a moved dependency closure
    Given a verified registered RefSeq snapshot with a selected genome FASTA
    When I prepare its exact versioned reference, snapshot ID and source path
    Then the original native bytes and paths remain unchanged
    And ordinary FAI and sequence dictionary files describe the complete FASTA
    And recipe identity is distinct from actual output SHA-256 values
    When I repeat the identical request
    Then verified prepared output is reused
    When I move the referenced snapshot and prepared directory together
    Then an independent indexed reader returns exact boundary and interior subsequences
    And complete records match an independent sequential reader
    And no BioV runtime, database, original path or network is required

  Scenario: Invalidate changed inputs and reject corruption
    Given a prepared index for one immutable snapshot
    When source content, implementation contract, tool version or parameters change
    Then the recipe identity changes
    When saved input or output bytes disagree with their recorded identities
    Then preparation reports corruption instead of returning a ready result
    And no existing native or prepared artifact is overwritten

  Scenario: Incomplete and duplicate work remain safe
    Given an interrupted preparation with partial staged output
    When a fresh preparation starts
    Then the staged output is not reused as a completed index
    And a verified complete index can be published without replacing a prior result
    When two identical preparations publish concurrently
    Then the verified immutable winner is retained

  Scenario: Reject unsupported FASTA without changing its meaning
    Given duplicate names, invalid wrapping or unsupported sequence bytes
    When preparation validates and indexes the input
    Then the request fails without publishing a prepared result
    And no input normalization, renaming or whole-sequence allocation occurs
```

D15: Native sequence-aware tables connect one already verified D14 RefSeq genome
FASTA preparation to D12 Rust Polars. The native `dataset_fasta_windows` tool
requires canonical exact reference, snapshot ID, selected native source path,
exact current recipe ID, exact opaque sequence ID and positive window width
(at most 1 MiB). Preview defaults to 5 rows and is bounded to 50. It never creates
a missing preparation, selects a latest snapshot/first sequence, downloads,
normalizes or changes native/prepared files. Established noodles indexed queries
extract bounded windows; native/prepared identities are verified before and after
reading, and conservative row/cell/identifier charges precede table allocation.
The independent 64 MiB retained-session charge and 16-dataset limit remain intact;
no process-wide peak-memory, large-data execution or speed guarantee is implied.

Complete non-overlapping windows use source-sequence-relative zero-based half-open
coordinates; a final partial row is retained. Ordered columns are `sequence_id`
(string), `start`, `end`, `length` (int64 bases), `is_full_window` (boolean),
`canonical_base_count`, `gc_base_count` (int64 bases), nullable `gc_fraction` and
`weighted_gc_fraction` (dimensionless float64). Canonical GC counts literal G/C
only and divides by case-insensitive A/C/G/T; all ambiguity is excluded and a zero
denominator is null. Weighted GC preserves existing core equal-base-set IUPAC
semantics, with every base in the denominator. Whole-sequence counts/fractions
are accumulated in the same pass and returned as `sequence_summary`; fractions
are not averages of row percentages.

The result is an ordinary dataset handle/schema/count/preview plus that summary.
Existing complete-data query/export/release and paired-record reopen work on the
typed table. Origin is structured lineage, separate from historical operation
handles: exact reference/snapshot/recipe, sequence ID/length, window width,
FAI/dictionary hashes, native FASTA hash/size, coordinate/GC policies and algorithm
revision. No CSV schema-policy claim is attached to FASTA. Its strict export/reopen record
is version 3, reserved for this validated sequence origin; existing CSV records
remain version 2 and record versions 1/2 remain readable. Portable companions
explain known column meanings without guessing species or provider versions.
A four-file window Arrow export is independent of raw source/prepared paths;
recomputation requires their dependency closure. Reopen preserves supplied origin
while verifying artifact consistency, never its authenticity or raw inputs.

### D15 acceptance cases

```gherkin
Feature: Exact native sequence metrics with complete portable table results
  Scenario: Scientific semantics and bounded previews
    Given an exact verified prepared mixed-case IUPAC DNA sequence
    When I request its positive-width windows by exact sequence and recipe IDs
    Then all windows use zero-based half-open source-relative coordinates
    And the final partial window is retained
    And all-ambiguous canonical GC is null while weighted GC is IUPAC-defined
    And the whole-sequence summary agrees with independent complete-base counts
    And a bounded preview is not substituted for complete table execution

  Scenario: Export, restart and independently move complete results
    Given a complete window dataset with more rows than the preview
    When I filter, sort and export through the existing dataset tools
    Then source identity, recipe, policies and operations remain discoverable
    When a new process reopens the saved record and Arrow pair without raw inputs
    Then a fresh handle preserves full typed rows and recorded lineage
    When the four-file export moves and historical directories disappear
    Then a standard PyArrow reader without BioV verifies full records and summaries

  Scenario: Exact selection and allocation failures remain read-only
    Given a missing/corrupt preparation or mismatched recipe/sequence selector
    When metrics are requested
    Then no new preparation or dataset handle is published
    Given a window request exceeding the remaining session row/cell memory charge
    When preflight rejects it
    Then existing data and all available handle slots are preserved
```

## §I INTERFACES
The following entries describe the existing Python interfaces unless marked
otherwise. The independent Rust MCP dataset route is specified in D12 and tracked
in T63–T67; the native dataset guide records its verified scope.

- native sequence tool: `dataset_fasta_windows` with configured native store → D15 exact prepared sequence windows, same-pass whole summary and D12 typed-table query/export/reopen; no Python analysis wrapper
- native prepared API: `biov_prepared::PreparedStore::new(store_root)` and `prepare_fasta(PrepareFastaRequest)` → D14; CLI `biov prepared fasta --store-root DIR --request-file JSON` and MCP `prepared_fasta` are thin bounded adapters
- native cmd: `biov storage register --store-root DIR --source-root DIR --request-file JSON` and `biov storage resolve --store-root DIR --request-file JSON` → thin adapters for D13, no downloads
- native tool: `storage_register` and `storage_resolve` through `biov mcp-native --data-root DIR --output-root DIR --store-root DIR` → the same bounded offline native-store contracts; `--store-root` is optional for the existing dataset-only route
- native library: `biov_storage::NativeStore::new(store_root)`, `register(source_root, RegisterRequest)` and `resolve(ResolveRequest)` → source-native immutable copy registration and offline filesystem resolution; D13 defines the bounded provider/validation scope
- file: `<store-root>/artifacts/<namespace>/<canonical accession>/snapshots/sha256-<digest>/` → relative inventory/receipt/checksums/README and complete unmodified `source/`; no durable index required
- api: `BioDataFrame.overlap(other, how, seqid_col, start_col, end_col, strand_col)` → selected self rows
- api: `BioDataFrame.intersect(other, seqid_col, start_col, end_col, strand_col)` → clipped self rows per overlap pair
- api: `BioDataFrame.subtract_ranges(other, seqid_col, start_col, end_col, strand_col)` → residual self fragments
- api: `BioDataFrame.nearest(other, seqid_col, start_col, end_col, strand_col, suffix, how)` → self + nearest other columns + `Distance`
- dtype: `biov.dna` | `biov.rna` | `biov.protein` → validated nullable uppercase sequence storage
- accessor: `Series.seq.length` → nullable integer Series
- accessor: DNA/RNA `Series.seq.reverse_complement()` | `gc_fraction()` | `translate(table=1, to_stop=False)`
- accessor: protein `Series.seq.molecular_weight()` | `isoelectric_point()` → nullable float Series; `amino_acid_composition()` → 20-column percentage DataFrame
- resource: `refseq.gcf://{+accession}` → original NCBI genome-summary JSON; `uniprot://{+accession}` → original complete UniProtKB JSON; other supported data namespaces → JSON file description (`uri`, `default_kind`, `representations`), without file download
- resource: `identifiers://{registry}` → native registry namespace object; `identifiers://{registry}:{+id}` → native resolver JSON
- api: `artifact_capabilities()` → independent JSON-compatible copy of supported namespaces, provider names, artifact kinds, defaults, and representation fields; no network or artifact-cache I/O
- file: `src/biov/assets/artifact_capabilities.json` → versioned single source for supported namespace/kind pairs, defaults, fixed file operations, and NCBI `fileType` selections
- tool: `parse_identifiers(prompt)` → JSON summary and ordered unique MCP resource links for resolver-valid explicit IDs, BioV resource URIs & allowlisted bare IDs
- tool: `resolve_identifiers(uri)` → the same resource content as a standard MCP embedded resource for tool-only clients
- cmd: `biov mcp` → BioV MCP server over stdio; ⊥ standalone `biov-mcp`
- api: `biov.analysis.run_analysis(AnalysisRequest)` → recorded synchronous Python analysis in a declared local Pixi environment; `inspect_analysis(record)` → bounded saved facts without resubmission
- tool: `run_analysis(request)` and `inspect_analysis(record)` → the same managed analysis path and structured summaries; failed runs retain MCP tool-error semantics, while querying a failed run can succeed
- resource: `file://{+path}` → complete registered analysis output, record or diagnostic bytes, subject to a 1 MiB resource-response cap and execution-root/identity checks
- cmd: `biov analyze REQUEST.json` and `biov inspect-analysis RECORD` → shared analysis execution and record inspection
- config: `BIOV_ANALYSIS_ROOT` → persistent outputs under the platform user data directory by default; optional `BIOV_ANALYSIS_BASE_URL` → existing HTTP(S) storage prefix for complete client downloads
- file: `docs/guides/analysis.md` → first two-step acceptance contract and verified client/platform scope; `docs/examples/sequence-analysis/` → ordinary scientific scripts and their native Pixi manifest/lock
- cmd: `biov update` → fetch & atomically refresh the packaged registry asset
- file: `src/biov/assets/identifiers_org_registry.json` → original complete resolver-dataset JSON response body
- api: `parse_identifier(value)` → exact single `IdentifierRef`; same generated URI, Compact ID, identifiers.org URL & curated bare-ID syntax as prompt parser
- api: `path(identifier, artifact=None)` → namespace-default environment-local immutable `Artifact` implementing `os.PathLike[str]`
- api: `open(identifier, artifact=None, mode="rb")` → handle for cached artifact
- fsspec: every manifest namespace has an installed read-only `BioVFileSystem` entry point; full identifier URIs retain namespace and select the same file as `path`; `artifact` is a filesystem/storage option
- api: `read_fasta(uri)` → existing dictionary of sequence records; `read_gff3(uri, storage_options={"artifact": "annotation_gff3"})` → existing `BioDataFrame` for a RefSeq annotation
- provider: `refseq.gcf` × `genome_fasta` → original catalog-selected genomic FASTA inside a complete NCBI Datasets package
- provider: `refseq.gcf` × `annotation_gff3` → original catalog-selected GFF3 inside a complete NCBI Datasets package
- provider: `refseq.gcf` × `rna_fasta` → original catalog-selected RNA FASTA inside a complete NCBI Datasets package
- provider: `refseq.gcf` × `cds_fasta` → original catalog-selected CDS FASTA inside a complete NCBI Datasets package
- provider: `refseq.gcf` × `protein_fasta` → original catalog-selected protein FASTA inside a complete NCBI Datasets package
- provider: `uniprot` × `protein_fasta` → original accession-named UniProtKB FASTA response
- provider: `uniprot` × `entry_json` → original complete accession-named UniProtKB JSON response
- provider: `uniprot` × `alphafold_cif` → one explicitly requested AlphaFold model
- provider: `pubmed|clinvar|dbsnp|geo` → available PMC PDF, complete VCV XML, complete RefSNP JSON, and single GEO matrix or explicit SOFT
- provider: manifest-declared structural/sequence/article/study/label/pathway namespaces → their native provider files; `encode` accepts only ENCFF file IDs
- api: `align_paired_reads(reference_fasta, fastq_r1, fastq_r2, bam_path, *, mem_args=(), sort_args=(), build_index=True)` → local BWA index/MEM execution through BioV software environments and pysam-sorted BAM; optional BAM index; native arguments and failures retained
- api: `crisprprimer` → existing deterministic matching, scoring, rice identifier conversion and fixed computation; BioV distributes the package and its lookup tables, while GEEPilot owns task selection and interpretation
- cmd: `crisprprimer`, `crisprprimer-docker`, `biov-azimuth` → migrated native Python workflow, existing Docker report bridge, and existing Azimuth backend runner; original scientific settings and scoring behavior retained
- file: `skills/` → specialized analysis and data-query guidance; no generic agent runtime
- env: uv + uv.lock → BioV development and Python 3.12/3.13/3.14 compatibility tests; scientific dependencies belong to deployment-selected Pixi environments
- cmd: `biov install [--environment-root DIR] [--pixi FILE] [--uv FILE] samtools|goatools` → Linux x86_64 thin upstream installation: Pixi 0.81.0 global installs Samtools 1.24 from conda-forge/bioconda and exposes its native command; uv tool installs GOATOOLS 1.6.5 plus statsmodels 0.14.6 using the complete BioV wheel's paired installed Python without downloading Python; managers resolve remaining dependencies rather than consume the bundled scientific lock; dedicated-bin PATH instructions are printed, never applied to the parent shell
- cmd: `biov list [--environment-root DIR] [--pixi FILE] [--uv FILE]` → official human-readable Pixi-global/uv-tool list output under path-labeled sections for existing dedicated roots only, failing on unreadable/non-UTF-8/over-1-MiB manager stdout rather than truncating; no custom JSON schema, parsed uv-text inventory or independent readiness/integrity audit; catalog/cache/locked execution-only environments and ordinary external manager installations are excluded
- cmd: `biov uninstall [--environment-root DIR] [--pixi FILE] [--uv FILE] samtools|goatools` → delegate removal of the selected owned backend tool environment and exposed commands; retain the other tool, caches, scientific data, outputs and separate locked workflow environments; no general purge/GC interface
- file: `environment_root/pixi-global/{envs,bin,cache,manifests/pixi-global.toml}` and `environment_root/uv-tools/{tools,bin,cache,python}` → isolated upstream-owned environments/entrypoints/native metadata; Pixi's owned global manifest is explicitly selected without user/XDG fallback; no BioV installed registry/journal, copied runner, custom launcher, bin override or hash-bound global version map; upstream commands survive removal of the management BioV wheel if their environments and GOATOOLS base Python remain available
- cmd: `biov python COMMAND ...` → explicit compatibility access to existing Python CLI capabilities, including manager-only and broader environment setup; native lifecycle routes remain Rust-owned
- cmd: `biov tools inspect|exec` → initial local Rust migration for bundled samtools/goatools on Linux x86_64, existing pinned Pixi 0.81.0, immutable content-keyed manifest/lock and successful native setup receipt; `exec --no-install` requires that receipt and a consistent Pixi prefix marker, passes literal native argv/stdio/status, and never installs; default exec performs locked setup or locked manager repair of an absent selected executable if needed; native execution uses validated typed Pixi activation and direct OS argv, including zero arguments; this is not the complete D6–D8 lifecycle, a project-manifest/SSH route or an MCP deployment surface
- cmd: `biov exec [--no-install] [--cwd DIR] [SOURCE:]NAME [ARGS]...` → bare NAME and conda:NAME select a declared environment's same-name Pixi task or executable, otherwise temporary Pixi execution; pypi uses uv tool run, npm uses npx; arguments pass through without requiring `--`
- config: `--config` or `BIOV_CONFIG` selects `config.toml` in the user configuration directory → application settings; environment variables override TOML; execution host uses native SSH configuration
- env: `BIOV_EXECUTION_HOST`, `BIOV_SSH_CONFIG`, `BIOV_EXECUTION_CWD` → private SSH destination/configuration and execution paths; no public remote subcommand
- env: `BIOV_ENVIRONMENT_MANIFEST` → explicit project manifest; unset selects bundled pyproject.toml and pixi.lock
- env: `BIOV_ENVIRONMENT_ROOT` → managed Pixi and scientific environments under the platform user data directory; used only when no configured or PATH Pixi matches the pinned version
- env: `BIOV_PIXI_BIN` → explicit Pixi executable; otherwise a matching `pixi` on PATH, otherwise the managed copy
- env: `BIOV_MAX_FILE_BYTES` → optional ceiling for each downloaded or decompressed provider file
- cmd: `biov python setup [ENVIRONMENT | --all] [--archive FILE] [--update-lock]` → make the pinned Pixi available (reusing a compatible one) and install a declared Pixi environment from the bundled manifest and lock
- cmd: `biov pixi [ARGS]...` → run the resolved Pixi manager with native arguments and its exit status
- cmd: `biov run [OPTIONS] SCRIPT [ARGS]...` → run complete Python script locally or submit it to LSF
- env: `BIOV_LSF_PYTHON` ? Python executable visible from LSF execution hosts; default = submitting interpreter
- file: `$BIOV_HOME/artifacts/refseq.gcf/<requested_accession>/` → unmodified extracted NCBI Datasets package root (`README.md`, `md5sum.txt`, `ncbi_dataset/...`)
- file: `$BIOV_HOME/artifacts/uniprot/<accession>/` → independently cached unmodified `<accession>.fasta` and/or complete `<accession>.json`; ⊥ manifest
- file: `.github/workflows/ci.yml` → test matrix, hooks & wheel/sdist inspection on push and pull requests
- file: Python sdist → package sources and resources, build metadata and licenses, Python tests, documentation and documentation build files; plugin sources are distributed through Git
- file: `.github/workflows/registry-drift.yml` → scheduled asset sync; opens a pull request when upstream changes

## §R RESEARCH
id|topic|finding|src
R1|RuRanges surface|stateless `ruranges.numpy` functions accept/return NumPy arrays; groups integer-coded; strand boolean only where required|https://github.com/pyranges/ruranges_py
R2|RuRanges kernels|`overlaps`, `nearest`, `subtract` return source indices; nearest physical directions = `forward`/`backward` & overlap distance = 0|https://raw.githubusercontent.com/pyranges/ruranges_py/master/ruranges/numpy.py
R3|pandas extension|custom dtype + 1-D ExtensionArray preserve semantic type; accessor init ! reject wrong dtype with `AttributeError`|https://pandas.pydata.org/docs/development/extending.html
R4|Biopython sequence math|weighted GC defines ambiguous IUPAC handling & empty GC = 0; molecular weight requires unambiguous residues|https://biopython.org/docs/latest/api/Bio.SeqUtils.html
R5|identifiers.org resolution|Compact Identifier = `[provider/]namespace:accession`; resolver `GET https://resolver.api.identifiers.org/{COMPACT_ID}` returns provider resources|https://docs.identifiers.org/pages/api.html
R6|MCP Python SDK|v2 `MCPServer` exposes typed tools & URI-template resources over stdio|https://py.sdk.modelcontextprotocol.io/
R7|identifiers.org registry dataset|`GET /resolutionApi/getResolverDataset` returns complete registry including namespace patterns, provider resources & institutions|https://docs.identifiers.org/pages/api.html
R8|NCBI genome CLI|`datasets download genome accession <GCF> --include ... --filename <zip> --no-progressbar` downloads one official ZIP; include values cover genome, GFF3, RNA, CDS, protein & sequence report|https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/
R9|NCBI genome package|extracted package root contains `README.md`, `md5sum.txt`, `ncbi_dataset/data/...`; each assembly keeps original files under its accession directory & catalog labels genome FASTA `GENOMIC_NUCLEOTIDE_FASTA`|https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/
R10|LSF submission|`bsub` acceptance assigns job ID while job may remain `PEND`; `DONE|EXIT` occur later ∴ ordinary submission ≠ completion|https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=management-job-lifecycle
R11|UniProt individual entry|`GET https://rest.uniprot.org/uniprotkb/<accession>.fasta` is the documented direct retrieval form for one UniProtKB entry in FASTA|https://www.uniprot.org/help/api_retrieve_entries
R12|UniProt structure links|one UniProtKB entry can expose many PDB cross-references with distinct methods, resolutions & chain coverage ∴ UniProt accession ≠ unique PDB coordinate file|https://rest.uniprot.org/uniprotkb/P42212.json?fields=xref_pdb
R13|UniProt complete entry|individual-entry REST retrieval supports accession-qualified `.json`; unfiltered response retains full entry metadata including database cross-references|https://www.uniprot.org/help/api_retrieve_entries
R14|NCBI genome summary|`datasets summary genome accession <GCF>` returns assembled-genome metadata as JSON without downloading a genome data package|https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/summary/genome/datasets_summary_genome_accession/
R15|PMC PDF access|the ID Converter identifies the current PMC version; cloud metadata supplies available PDF object and MD5; old OA web service was retired in August 2026|https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/
R16|GEO files|Series matrices contain values with series/sample metadata; full SOFT preserves accession data; multiple matrices need explicit selection|https://www.ncbi.nlm.nih.gov/geo/info/download.html
R17|ClinVar complete record|EFetch with rettype=vcv and is_variationid retrieves the complete variation XML|https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/
R18|dbSNP full data|RefSNP endpoint returns the full native variation JSON by rs number|https://api.ncbi.nlm.nih.gov/variation/v0/

R19|uv tool lifecycle|uvx is uv tool run; cached run environments are disposable, installed environments persist until uninstall; explicit versions, refresh and upgrades govern reuse|https://docs.astral.sh/uv/concepts/tools/
R20|uv storage|persistent tools, disposable cache and command directories are distinct and inspectable/configurable|https://docs.astral.sh/uv/reference/storage/

## §V INVARIANTS
V1: intervals use 0-based, end-exclusive `[start,end)`; integer `0 ≤ start < end`; touching boundaries ≠ overlap
V2: operations group by exact `seqid`; when `strand_col` exists on both, only `+|-` valid & group by exact `seqid,strand`; one-sided/invalid strand → `ValueError`; `strand_col=None` ignores strand
V3: outputs preserve self input order; pair expansions use self row then other row order; fragments use self row then ascending coordinate; duplicate input rows remain distinct
V4: empty self/other inputs return typed empty or unchanged results with stable schema; ⊥ kernel panic
V5: `overlap` returns each matching self row once; `first|last` select membership, `containment` means self contains other, `member` means self contained by other
V6: `intersect` emits one clipped self-metadata row per overlap pair; overlap duplicates yield duplicate output rows
V7: `subtract_ranges` emits every non-empty residual fragment with self metadata; overlapping/duplicate masks do not duplicate residual space
V8: `nearest` returns ≤1 row/query; overlap `Distance=0`, otherwise `gap+1` so adjacent half-open intervals have `Distance=1`; `next|previous` = genomic right|left; `upstream|downstream` = strand-aware 5′|3′; tie → lowest other input row; missing group candidate → omit query
V9: BioV core dependencies exclude PyRanges, sorted-nearest and ncls; optional analysis environments follow their scientific packages' native dependencies
V10: ordinary string Series `.seq` → `AttributeError`; caller ! choose/carry `biov.dna|rna|protein`; ⊥ content guessing
V11: sequence storage uppercases valid strings, preserves `pd.NA`, accepts empty strings & declared IUPAC alphabets, rejects non-string/invalid symbols with `SequenceValidationError`
V12: DNA/RNA reverse complement preserves dtype/nulls; weighted GC handles IUPAC & empty string; translation requires complete codons, honors `table,to_stop`, preserves nulls, returns `biov.protein`
V13: protein mass/pI/composition use non-empty canonical 20 amino acids; extended IUPAC, stop, or empty sequence stored but analysis → stable `SequenceValidationError`; null results remain null; composition columns = canonical amino-acid order & values = percentages
V14: the existing Python implementation confines RuRanges calls to the interval module and remaining Biopython sequence algorithms to the sequence module; migrated normalization/reverse complement/length/weighted GC use the Rust sequence core; D9–D10 supersede backend exclusivity for the target implementation, which centralizes migrated rules in Rust under the retained or explicitly revised contracts
V15: preserve existing behavior except explicitly revised contracts; artifact access uses standard exceptions without BioV-specific error subclasses
V16: prompt parser recognizes explicit Compact Identifiers, identifiers.org URLs, BioV resource URIs & allowlisted bare IDs; canonicalizes variants, preserves first occurrence order, deduplicates, accepts provider/slash accessions & trims prose punctuation
V17: parser emits links only for resolver-valid IDs; invalid candidates omitted; upstream/service/JSON failures surface distinctly ≠ invalid ID
V18: `identifiers://registry:id` round-trips reserved accession characters & returns native resolver JSON; `identifiers://registry` losslessly reconstructs the native namespace object; invalid registry/accession → stable resource error
V19: identifiers.org integration is read-only/idempotent, accesses only the fixed resolver, registry-dataset and identifiers.org origins, has bounded prompt candidates & request timeout
V20: the existing Python `biov mcp` starts stdio without protocol-corrupting stdout; `biov-mcp` ∉ installed scripts; expected provider errors are translated once to SDK errors so clients receive their messages
V21: asset preserves complete official response in native nested shape; namespace resources, institutions & locations remain upstream-owned objects
V22: MCP publishes a data template for each manifest namespace, the two `identifiers` forms; no API catalog or old `identifiers://resolve/*` resource
V23: data template metadata exposes prefix, accession regex, sample & namespace ID; every ID read validates the asset regex; `namespaceEmbeddedInLui` reconstructs resolver input correctly
V24: prompt tool uses resolver-parsed namespace/local ID → direct resource for a manifest namespace, otherwise generic identifiers resource; provider-qualified input maps to a location-independent URI; output URI deduplicated in prompt order
V25: asset validation checks only fields required for runtime indexing/routing; unknown upstream fields & nesting ! preserved; asset ! packaged in wheel/sdist
V26: registry update validates upstream JSON/runtime fields; byte-identical body skips replacement; changed body atomically replaces asset unchanged
V27: bare-ID recognition checks only an explicit namespace-prefix allowlist; `refseq.gcf` ∈ allowlist; non-allowlisted registry patterns ∉ inference; invalid shorthand omitted before resolver access
V28: each accepted identifier syntax variant ! pass an isolated acceptance case without another valid fallback form; combined variants ! canonicalize & deduplicate before resolver access
V29: `parse_identifier` accepts exactly one full reference, validates packaged namespace regex, returns direct data URI for supported namespaces or generic identifiers URI otherwise without resolver I/O; prose, multiple references, unsupported per-registry schemes, unknown namespaces & non-allowlisted bare IDs → stable syntax error
V31: versioned GCF path request returns that exact catalog assembly; versionless GCF accepts exactly one matching versioned catalog assembly; missing/ambiguous/non-GCF package → FileNotFoundError
V32: `Artifact` is accepted wherever `os.PathLike[str]` is accepted & exposes original package member path, package root, requested/canonical identifier, kind & byte size
V33: valid package-directory cache hit reads only local catalog/file metadata & skips `datasets`; cache miss downloads/extracts in a unique sibling staging directory then atomically publishes the complete package root; failed/interrupted writes never become cache hits
V34: command argv = `datasets download genome accession <accession> --include gff3,rna,cds,protein,genome,seq-report --filename <temporary-zip> --no-progressbar`; every safe ZIP member retains its exact relative path; ⊥ renamed/copied FASTA, BioV manifest, partial extraction or direct REST download
V35: original package catalog selects exactly one existing `GENOMIC_NUCLEOTIDE_FASTA`; unsafe/duplicate ZIP paths and missing/duplicate FASTA → ValueError; accession mismatch → FileNotFoundError; native parser and I/O errors propagate; rejected `datasets` CLI reports its diagnostic
V36: `biov run` passes Python executable, absolute script & arguments as argv without shell interpolation; complete script—including `path`—runs inside selected executor
V37: local executor inherits current cwd/environment/stdio & command exit code becomes `biov run` exit code
V38: LSF executor submits via `bsub`, pins cwd, supports queue/name/stdout/stderr + executor-visible Python, returns parsed numeric job-ID receipt; accepted submission never claims job completion; missing/rejected/unknown receipt → stable execution error
V39: cached path is valid only inside current executor; LSF submission does not resolve or return compute-node paths to submitter; shared env/cache/network availability remains deployment configuration
V40: unsupported namespace × artifact, invalid identifier, assembly mismatch, malformed package, missing/rejected `datasets`, missing executor & rejected submission have distinct public exception types/messages
V41: omitted artifact kind dispatches by the manifest namespace default; explicit unsupported pairs raise `ValueError` before download
V42: UniProt request URLs = fixed HTTPS origin + `/uniprotkb/<percent-encoded-accession>.fasta|.json`; streamed response bytes remain unchanged at accession-named files; ⊥ generated manifest, sequence/JSON rewrite or field filtering
V43: UniProt FASTA/JSON validates its local regular file, bounded FASTA header or JSON object, and exact accession; a miss atomically publishes only the requested representation; mismatched content never becomes a cache hit
V44: `uniprot://<accession>` defaults to protein sequence; PDB coordinates require an explicit PDB ID; explicit `alphafold_cif` requires exactly one matching model, with no arbitrary selection
V45: `entry_json` and `protein_fasta` cache independently; requesting one never fetches or requires the other; full raw JSON retains all upstream PDB IDs and cross-reference properties
V46: registry asset bytes = upstream response body; ⊥ wrapper, flattening, foreign keys or derived fields; runtime indexes native records in memory without mutating them
V47: `refseq.gcf://<accession>` executes only `datasets summary genome accession <accession>`; validates one matching report then returns stdout unchanged; ⊥ `path`, package download/extraction or `$BIOV_HOME` write; `path(refseq.gcf)` remains V31–V35
V48: each registered `refseq.gcf` artifact kind maps to exactly one official catalog `fileType`; the catalog selects exactly one existing member per request; missing/duplicate members → ValueError; native parser and I/O errors propagate; one cached package serves every kind
V49: provider integration publishes no API catalog assets, catalog resources, or generic upstream query dispatcher; search and scientific API queries use task skills and upstream OpenAPI/GraphQL documentation, never pickle schemas
V50: file namespaces validate accessions with the packaged identifiers.org regex before provider I/O; syntax validity does not imply file availability
V51: file operations are explicit manifest representations or bounded provider lookups; callers cannot supply arbitrary HTTP endpoints through the artifact API
V52: PubMed uses the current PMC version and checks its PMID, version, PDF object, and published MD5; no available PDF → failure, never an abstract or scraped-page substitute
V53: ClinVar downloads complete matching VCV XML; dbSNP downloads matching complete RefSNP JSON; GEO defaults to `expression_matrix`, requiring exactly one nonempty matching GSE Series matrix and rejecting missing or multiple matches; native full SOFT is available only when explicitly selected
V54: a native file format does not imply completed biological normalization or analysis; multiple GEO matrices, ENCODE experiments, and multiple AlphaFold models require explicit caller decisions
V55: existing Python MCP file descriptions return URI, default kind, and representations without file download or cache I/O; RefSeq/UniProt retain V47/V45 metadata behavior; actual files are accessed through Python or fsspec
V56: the versioned artifact manifest owns namespaces, kinds, defaults, and NCBI selection fields; discovery validates dispatch structure without network/cache access; caller mutations cannot change dispatch; installed fsspec schemes match manifest namespaces
V57: individual downloads stage and validate before atomic publication; valid local cache reads skip downloading; format, identifier, HTTP, and network failures do not become successful error-text results or fallback files
V58: fsspec discovers every manifest namespace without prior `import biov`; reads and `info` reuse `path` validation/download/cache; handles support seek and ranged reads, `open_local` returns the actual filename; writes fail before provider I/O, errors propagate, and directory listing remains unsupported
V60: local cache roots come from BIOV_HOME, native fsspec configuration, or explicit directory arguments; BioV preserves existing fsspec filecache settings; packaged metadata and repository examples contain no workstation-specific cache paths
V62: `biov exec` treats only conda, pypi and npm prefixes as sources; bare NAME defaults to conda; bare NAME and conda:NAME select a declared Pixi environment of the same name if present, otherwise `pixi exec -s NAME -- NAME ARGS`; pypi:NAME delegates to `uv tool run NAME ARGS`, npm:NAME to `npx --yes NAME ARGS`; declared environments use `install --locked`, their declared preparation task once, then `run --frozen`, resolving a same-name Pixi task before a same-name executable; BioV options are parsed only before the coordinate, and native arguments reach the program unchanged after an optional initial `--` separator, including quotes, empty strings and shell metacharacters; `--no-install` uses `run --as-is` for declared environments and rejects temporary sources before execution; missing prerequisites fail with actionable errors without switching sources or falling back to host PATH; stdio and native status are preserved
V63: commands from every source use execution-host paths and inherit caller caches; BioV supplies BIOV_HOME and BIOV_CACHE_HTTP; package/environment configuration uses native manager controls, with no per-software configuration table or cache rewriting
V64: runtime/environment configuration belongs to deployment; skill instructions expose scientific programs and native arguments; MCP adds no deployment-management tools; dependency checks do not imply scientific validation
V65: SSH is internal to exec and delegates authentication/configuration/host-key checking to OpenSSH; argv is POSIX-shell-quoted; only the coordinate, native arguments, optional cwd and no-install selection are forwarded; remote manager/cache settings belong to the remote configuration; no automatic file transfer or public remote subcommand
V66: setup stays explicit and local to the execution host; a configured, PATH or managed Pixi is reused only when its reported version matches the pinned release, and only a downloaded archive must match platform SHA-256 and retain the Pixi license; scientific dependencies and versions belong to src/biov/assets/environments/pyproject.toml and its adjacent pixi.lock, shipped in the wheel and copied into writable content-addressed workspaces when no project is selected; setup uses install --locked and rejects missing/stale locks; exec provisions declared environments from that same lock on demand and never rewrites it, while temporary-source execution resolves packages through its native manager without that project lock; `--no-install` skips provisioning for declared environments and is rejected for temporary sources; shipped channels exclude Anaconda defaults; machine settings remain in application TOML and environment/package declarations use native manager files; preparation tasks an environment declares for source-only tools run on setup and once on demand after installation; initialization requires no agent or LLM
V67: the sdist contains the Python package and required resources, build metadata and licenses, Python tests, documentation and its build files; source selection is declared in MANIFEST.in and pyproject.toml and CI checks the selected project-file set exactly; skills, plugin manifests and plugin-only tests remain in the Git distribution; an unpacked sdist can build the wheel, run its Python tests and build the documentation
V68: declared environment entry points belong to native Pixi same-name tasks when the executable name differs; shipped entry tasks preserve the caller directory through INIT_CWD and locate prepared sources through PIXI_PROJECT_ROOT; library entries execute Python scripts, not invented library CLIs; a missing same-name task and executable produces an actionable error naming the environment and manifest, while an existing entry's own exit 127 is preserved; --no-install on an uninitialized bundled workspace reports biov python setup NAME without creating it; tasks and executable discovery use Pixi's resolved task and activation data
V69: managed analysis executes an ordinary Python script synchronously through existing `biov exec` in a declared local locked Pixi environment whose entry accepts scripts; configured SSH managed execution fails explicitly; scientific algorithms and task-specific orchestration remain in caller/example scripts
V70: every managed run has a unique private persistent directory with exact code, parameters, input copies, complete outputs and logs; inputs are copied with change checks, made read-only, hashed and checked after execution; prior-result reuse verifies the completed record and output identity; missing/modified results fail without recomputation, and failed runs never publish partial files as completed outputs
V71: the driver alone writes its atomic run record and the worker separately records confirmed scientific-process facts; launch intent precedes execution, launcher creation is not scientific startup, and unconfirmed post-launch failures or nonterminal records queried after restart report unknown without PID-based guesses or resubmission; program completion, declared checks and output validation remain distinct; second-step failure does not alter the first run's saved results
V72: default analysis summaries include method name, inputs, parameters, requirements, performed checks and at most two preview rows/records, with a 32 KiB total response cap and explicit omissions; complete files remain on disk, small registered files are available as bounded MCP resources, and a configured existing HTTP(S) storage endpoint supplies full-download URLs without automatic publication or a custom transfer protocol; supported scope and actual client verification are recorded in the analysis guide

V73: BioV owns the migrated deterministic crisprprimer and Azimuth implementations and verified publication assets; public Python imports and CLI behavior remain available without GEEPilot; the unpublished NAU mapping and its unused globals are excluded pending source/provenance verification, while public RAP/MSU conversion remains available; unsupported RAP annotation lookup fails explicitly; BLAT uses run_software with literal native argv and propagates failure before reading outputs; caller-local temporary files are never sent to configured SSH execution; cache defaults use BIOV_HOME or native fsspec configuration; scientific methods, model provenance, safety routing and interpretation remain task-skill responsibilities, and migration checks do not imply biological validation

V74: paired-read alignment uses run_software for BWA indexing and alignment, keeps local reference/FASTQ/SAM/BAM paths on the same host, retains caller-specified alignment and sort arguments, and propagates failed commands before consuming their outputs; fixture checks execute actual BAM sorting/indexing without claiming a live BWA run

V75: Rust migration defines and tests its selected CLI/MCP and any optional Python contracts; Rust-side Polars/MCP is the primary analytical direction and Arrow IPC provides table interoperability; breaking API/type changes are explicit and need no Python mirror or legacy shim; complete scientific outputs, metadata and data integrity remain validated independently of interface compatibility
V76: native release support is claimed only for tested distribution targets; no compiler is required for supported prebuilt wheels, while source builds declare their toolchain requirements
V77: every cache/data design obeys D0; copied/moved bundles must be discoverable, interpretable and analyzable with ordinary readers without BioV; preserve complete native data, explicit unknowns, relative inventory, content identities and lineage; current legacy gaps stay visible and require separate migration acceptance

V78: native copy-registration preserves complete source bytes/layouts and original ownership, publishes only validated immutable snapshots without replacement, and resolves offline through portable filesystem records; exact biological references, snapshot content identity, declared scope and representation availability remain separate; full checksums do not authenticate provenance, and stream-copy storage does not expand D12 analytical limits

V79: prepared reference indices preserve native bytes, use established format implementations, distinguish recipe identity from output hashes, verify reuse and no-replace publication, and declare relative dependency closure for independent moved-bundle readers; incomplete staging, corrupt input/output and unsupported FASTA never become ready artifacts

V80: user-tool lifecycle adapters use only explicit BioV-root-owned upstream state: Pixi 0.81.0 global Samtools 1.24 and uv tool GOATOOLS 1.6.5 plus statsmodels 0.14.6 with paired installed Python; native manifests/receipts and human-readable list output remain authoritative; no custom install registry, runner, journal, JSON list, bin override or full transitive lock is introduced; no user-global manager setting or shell profile is changed; removal targets only the selected backend tool environment/commands and preserves caches, scientific data/results and separate locked workflows; base Python availability is required after management-wheel removal

## §T TASKS
id|status|task|cites
T1|x|write contracts & failing acceptance tests|I.*,V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13
T2|x|replace interval adapter with RuRanges NumPy kernels|I.overlap,I.intersect,I.subtract_ranges,I.nearest,V1,V2,V3,V4,V5,V6,V7,V8,V9,V14
T3|x|add typed sequence EA/dtypes/accessor|I.dtype,I.accessor,V10,V11,V12,V13,V14
T4|x|update public docs, exports, dependencies & lock|V9,V15
T5|x|run full verification matrix|V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15
T6|x|write identifiers.org MCP acceptance tests|I.resource,I.tool,I.cmd,V16,V17,V18,V19,V20
T7|x|implement resolver client, resource, parser tool & stdio server|I.resource,I.tool,I.cmd,V16,V17,V18,V19,V20
T8|x|document MCP setup, parsing scope & error behavior|I.resource,I.tool,I.cmd,V16,V17,V18,V19,V20
T9|x|run full verification matrix|V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20
T10|x|write raw registry asset & generated-resource acceptance tests|I.resource,I.file,V18,V21,V22,V23,V24,V25
T11|x|persist raw registry response & generate full asset|I.file,V21,V22,V25
T12|x|generate namespace MCP resources & prompt links from asset; remove old URI|I.resource,I.tool,V16,V17,V18,V19,V22,V23,V24
T13|x|document registry URI schemes, aliases & asset refresh|I.resource,I.file,V22,V23,V24,V25
T14|x|run full verification matrix & package inspection|V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20,V21,V22,V23,V24,V25
T15|x|add `biov update` command & acceptance tests|I.cmd,I.file,V21,V25,V26
T16|x|replace standalone `biov-mcp` with `biov mcp` & update docs/package|I.cmd,V15,V20
T17|x|support generated URI variants & allowlisted unambiguous bare IDs in `parse_identifiers`|I.tool,V16,V17,V19,V24,V27,V28
T18|x|write single-ID, path resolution/cache/package & local/LSF execution acceptance tests|I.api,I.provider,I.cmd,V29,V31,V32,V33,V34,V35,V36,V37,V38,V39,V40
T19|x|implement `IdentifierRef`, artifact provider registry, initial NCBI GCF FASTA provider/cache & public exports|I.api,I.provider,I.file,V29,V31,V32,V33,V34,V35,V39,V40
T20|x|implement whole-script local/LSF executors & `biov run`|I.cmd,I.env,V36,V37,V38,V39,V40
T21|x|document LLM/Biopython usage, executor boundary, artifact cache & deployment requirements|I.api,I.cmd,I.env,V31,V32,V38,V39,V40
T22|x|run full tests, hooks, package inspection & bounded official NCBI smoke validation|V15,V29,V31,V32,V33,V34,V35,V36,V37,V38,V39,V40
T23|x|replace REST/flattened-FASTA acceptance tests with official CLI argv, complete-package layout, catalog selection & cache tests|I.api,I.provider,I.file,V31,V32,V33,V34,V35
T24|x|replace RefSeq REST adapter/manifest cache with `datasets` CLI & verbatim package cache|I.api,I.provider,I.file,V31,V32,V33,V34,V35
T25|x|document CLI prerequisite/package layout & run full verification|V15,V31,V32,V33,V34,V35
T26|x|write UniProt URI, official endpoint, raw-file cache & failure acceptance tests|I.api,I.provider,I.file,V27,V41,V42,V43,V44
T27|x|implement namespace-default artifact dispatch & UniProt FASTA provider|I.api,I.provider,I.file,V41,V42,V43,V44
T28|x|document UniProt/PDB boundary, run full verification & download official GFP `P42212`|V15,V41,V42,V43,V44
T29|x|backprop FASTA-only cache bug; preserve & expose complete UniProt JSON beside FASTA|I.api,I.provider,I.file,V42,V43,V44,V45
T30|x|replace relational registry snapshot with raw upstream response|I.file,V21,V25,V26,V46
T31|x|decouple RefSeq MCP summary reads from complete analysis-package path resolution|I.resource,V23,V47
T32|x|harden review findings: default MCP resource security, mapped registry CLI errors, bounded datasets/bsub subprocesses|V19,V20,V35,V38,V40
T33|x|generalize RefSeq catalog selection to annotation/rna/cds/protein artifact kinds|I.provider,V31,V32,V33,V35,V48
T34|x|add CI test/hook/package workflows & scheduled registry-drift sync; restore byte-exact asset|V15,V25,V26,V46
T35|superseded|initial API catalogs and NCBI summaries replaced by identifier file providers and task skills in T39|V49,V52,V53,V55
T36|x|implement manifest-driven artifact discovery and selection; synchronize migration inventory and validate package contents|V48,V56
T37|x|register read-only identifier filesystems with fsspec; preserve existing FASTA/GFF reader results and verify cache reuse and installed entry points|I.fsspec,V15,V58
T39|x|replace API catalogs and summaries with identifier file providers; retain task guidance in skills and synchronize documentation|V49,V50,V51,V52,V53,V54,V55,V56,V57,V58
T40|x|verify combined migration tests, package contents, dynamic fsspec entry points, and preserved existing behavior|V15,V25,V56,V58
T41|x|preserve configurable local caches and native fsspec settings; verify environment-selected cache paths|V60
T43|x|add native scientific-program execution in existing environments and execution-host file access; verify argv, file output and failure status|I.cmd,I.env,I.file,V60,V62,V63,V64
T44|x|make exec the sole native-command entry point; read application and SSH configuration internally, support arbitrary commands in host/pixi environments, and verify routing, literal argv, caches and remote status without implicit file transfer|I.cmd,I.env,V60,V62,V63,V64,V65
T45|x|use native Pixi manifests and locks, use remote-local configuration, inherit caller caches, and explicitly prepare scientific environments with pinned Pixi and community sources|V60,V62,V63,V65,V66
T46|x|align exec with uvx/npx: resolve a declared environment by name, provision it on demand from the lock including its declared preparation task, and add `--no-install`|I.cmd,V62,V66
T47|x|replace software aliases and runtime defaults with package coordinates across conda, PyPI and npm; default bare names to conda and preserve declared locks, preparation, caches, literal arguments and remote status|I.cmd,I.env,V62,V63,V65,V66
T48|x|include Python sources and resources, Python tests, documentation and its build files in the sdist; verify the exact selected file set, wheel build, unpacked tests and strict documentation build|I.file,V67
T49|x|restore declared tool entry points as native Pixi tasks, preserve literal task arguments and caller directories, and diagnose missing entries without masking native failures|I.cmd,V62,V68
T50|x|repair the pinned DiffDock inference environment and default configuration path; verify DiffDock's CLI and scvi's actual package import in locked Linux environments|I.cmd,V64,V66,V68
T51|x|fix the complete-CDS translation → protein-properties example and acceptance cases in docs/guides/analysis.md, including scientific checks, concrete interfaces, reference context, previews, retrieval, input stability and interruption behavior|D1,D2,D3,D4,D5
T52|x|implement shared local managed Python analysis through existing Pixi execution, with persistent records, complete outputs, bounded previews, verified references, MCP/CLI access and explicit HTTP retrieval|D1,D2,D3,D4,D5,V69,V70,V71,V72
T53|x|validate the real five-protein two-step example in locked macOS ARM64 Pixi, SDK stdio and Codex CLI 0.153.4, plus scientific and SDK stdio/HTTP acceptance in Docker linux/amd64 emulation; verify preview-independent reuse, failed-step retention, missing/changed outputs, interrupted records and complete HTTP downloads above the resource limit; record the exact tested scope in docs/guides/analysis.md; desktop and other clients are not claimed as validated; managed SSH/LSF is not implemented, and real stdio verification confirms configured remote analysis is rejected before creating results|D2,D3,D4,D5,V69,V70,V71,V72

T54|x|migrate existing GEEPilot deterministic CRISPR computations, Docker report parsing, Azimuth runner and assets into the BioV distribution; replace removed BLAT API with native execution and verify imports, native argv, failure propagation, cache reuse, inputs and package contents|V60,V62,V65,V67,V73

T55|x|migrate shared GEEPilot BWA-to-BAM execution into a public BioV API, remove the external PATH requirement, and test native argv, sort/index results and failure handling with SAM fixtures|I.api,V62,V65,V74

T56|partial|single native biov entry and thin isolated upstream install/list/uninstall contract defined for Samtools/Pixi-global and GOATOOLS/uv-tool on Linux x86_64; preserve setuptools-rust mixed packaging, paired-interpreter Python bridge and separate bundled locked scientific workflows; explicit upgrades and general cache cleanup remain gaps|D6,D7,D8,V80
T57|partial|replace custom installed records/launchers/retained runners/journals with isolated native Pixi-global and uv-tool adapters, primary-version pins, official raw human inventory and backend-selected removal preserving scientific data and locked workspaces; validate real managers/tools, isolated roots, failures/collisions, management-package uninstall survival with retained base Python, package upgrade and wheel/sdist installation before claiming completion; no fully locked transitive global-install or general upgrade/cleanup promise|D2,D3,D4,D6,D7,D8,V80

T58|partial|sequence normalization/reverse-complement public types, errors, breaking Unicode correction and full-output independent fixtures fixed in docs/guides/sequence-contract.md; inventory remaining scientific and operational contracts per slice, without requiring a full parallel Python analysis API|D9,D10,V75
T59|partial|core/PyO3/development binary workspace, native normalization/reverse complement implemented; original extension wheels/source builds validated; migration to setuptools-rust mixed executable/extension packaging requires fresh wheel/sdist acceptance; portable Linux/macOS release-target validation pending|D9,D10,D11,V10,V11,V12,V67,V75,V76
T60|partial|native validated sequence lengths and weighted IUPAC GC implemented with independent base-set fixtures; move remaining sequence algorithms, interval operations and format parsing to Rust in contract-sized slices; use Rust-side Polars for analytical tables, evaluate scientific library candidates, measure interoperability costs, document schemas and metadata handling and independently verify any owned algorithms|D9,D10,D12,V1,V2,V3,V4,V5,V6,V7,V8,V12,V13,V75
T61|partial|curated offline GCF/UniProt identifiers, local dataset provenance and the bounded native samtools/goatools locked setup-to-exec bridge are implemented; move remaining offline identifier parsing/validation, provider resolution/data provenance, cache, configuration, tool lifecycle/execution and managed results into their Rust responsibility boundaries; instantiate crates only when actual responsibilities justify them; retain external tool/native lock ownership and all failure/receipt/result contracts; implement T56–T57 lifecycle scope there|D2,D3,D4,D6,D7,D8,D9,D10,V75
T62|partial|the local native dataset MCP slice and source-installed binary are validated; expand native CLI/MCP scope after selected operation and real-client acceptance; keep legacy Python and native Rust routes explicit until any switch is independently justified; publish validated binaries and optional wheels/source builds, retiring superseded Python logic/dependencies when no retained function needs them|D9,D10,D11,D12,V20,V64,V67,V75,V76
T63|x|fix the first local Rust/Polars MCP dataset contract, string-safe defaults and explicit schemas, supported operations, handle lifetime, path boundaries and resource limits; distinguish unknown biological metadata and session handles from identifiers; record complete-data/export/chunk-retrieval/readback acceptance and the curated offline identifier subset|D2,D3,D4,D9,D10,D12,V75
T64|x|implement official rmcp stdio behind biov mcp-native with configured data/output roots, local CSV open/schema/row count/bounded preview, complete-data filter/sort/select into derived datasets, Arrow IPC plus record export, explicit bounded retrieval and dataset release; keep the existing Python MCP route separate|D9,D11,D12,V20,V75
T65|x|validate the new native workflow through an MCP client and independent Arrow IPC readback, including records beyond previews, identifier/numeric-text fidelity, explicit numeric schemas, duplicates/nulls/order, digest-verified chunk reconstruction, unknown metadata, handle expiry, errors/limits/path confinement and clean stdio; update the actual supported scope only after checks pass|D2,D3,D4,D10,D11,D12,V75,V76
T66|x|implement dataset_reopen for a previously exported root-confined JSON record and same-directory IPC pair, validating record version/format, bounded same-snapshot bytes and digest, schema/actual types, row count and metadata before issuing a fresh session handle; expose and persist reopening verification context without authenticating recorded provenance or biological metadata|D2,D3,D4,D9,D12,V75
T67|x|validate complete cross-process reuse in source-built Linux through real MCP processes A and B: typed CSV open/filter/sort/select/export/EOF, reopen/preview/query/re-export and independent full-table/digest readback; cover types/nulls/order/duplicates, expired handles, missing/tampered/mismatched pairs, limits and path escapes with no fallback; 29 biov-data tests and all nine real MCP subprocess cases pass|D2,D3,D4,D10,D11,D12,V75,V76

T68|x|add native export README and versioned companion manifest without changing strict record v2; emit upstream-compatible LargeUtf8 for direct standard-reader filtering and retain prior Utf8View paired reopen through bounded metadata/allocation preflight; validate moved-directory independent standard-reader checksums/schema/full rows/filter/summary with BioV absent, explicit unknown meanings/versions and honest trust scope; installed Linux no-cast PyArrow 25.0.1 acceptance, 39 biov-data tests and all nine MCP subprocess cases pass|D0,D2,D3,D4,D12,V77
T69|planned|migrate audited Python provider/fsspec/managed-analysis/CRISPR cache and output gaps to D0 in contract-sized slices; preserve original provider bytes and native layouts, add missing portable inventory/dictionaries/source identities/lineage and independent-reader acceptance; do not claim global compliance from T68|D0,D2,D4,D7,D10,V77

T70|x|implement the D13 bounded Rust native-store library and validate exact RefSeq catalog/checksum registration, explicitly declared PDB scope, offline paths, missing/ambiguous/unavailable/corrupt outcomes, no-overwrite concurrency, staging isolation, moved-store discovery and ordinary-reader portability; source-built Linux acceptance passes with 32 storage tests, 17 MCP/CLI process cases and independently read relocated real RefSeq/PDB packages under kernel network blocking; see native-storage guide for limits|D0,D9,D13,V77,V78
T71|planned|add further native semantic adapters and analysis-ready derived views only through independent provider/format contracts; separately gate transactional downloads/hydration, import/result lineage, optional indices, explicit GC, large-data execution and distributed storage; five-source observations are not implementation claims|D0,D9,D10,D13,V77,V78

T72|x|implement and independently validate D14 pinned RefSeq genome FASTA preparation with upstream noodles-fasta, conventional FAI/TSV outputs, bounded streaming, recipe invalidation, verified reuse, atomic publication and moved offline standard-reader acceptance; source-built Linux passes 25 prepared-core tests, 21 installed CLI/MCP cases and 10 independent portability cases including actual RefSeq and SIGKILL/retry; see prepared guide for adversarial coverage and explicit platform/input limits|D0,D9,D13,D14,V77,V78,V79

T73|x|implement D15 exact prepared RefSeq sequence-window GC tables and same-pass whole counts, conservative allocation preflight, structured origin/known column meanings and existing full-query/export/reopen integration; all four installed Linux acceptance cases pass, including offline moved standalone PyArrow 25.0.1 analysis without BioV, actual E. coli complete 465-window counts and 400001-row rejection preserving all handle slots; see prepared guide for bounded scope and core validation|D0,D9,D10,D12,D13,D14,D15,V75,V77,V79

## §B BUGS
id|date|cause|fix
B1|2026-08-30|mixed-form acceptance prompt contained Compact ID fallback, masking absent URI/bare parsing|V28
B2|2026-08-30|new allowlist error path exposed incomplete MCP helper docstring contracts|docstrings completed; no behavior invariant
B3|2026-08-31|RefSeq provider reimplemented NCBI REST and flattened one renamed FASTA instead of retaining the official CLI data package|V31,V32,V33,V34,V35
B4|2026-09-01|UniProt provider conflated default computation artifact with complete upstream cache & downloaded FASTA only|V45
B5|2026-09-03|registry snapshot transformed upstream JSON into an unused relational format|V46
B6|2026-09-03|RefSeq MCP metadata read reused `genome_fasta` path resolution & downloaded complete analysis package|V47
B7|2026-10-01|T47 removed command aliases without moving differing executable and prepared-source entry points to native tasks; run --executable also bypassed tasks|V62,V68
B8|2026-10-01|T49 moved differing entry points to native tasks but missed `lumpy`, whose environment name matched a different real binary, so exec ran the low-level caller without failing|V68
B9|2026-10-01|exec diagnosed a 127 exit through Pixi introspection that aborted on error, so a failed diagnosis replaced the entry's own status instead of preserving it|V68
