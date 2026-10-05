# Selected model resources

`biov model download` retains an explicit selection of native files from a
Hugging Face **model** repository at a caller-supplied full Git commit. The
official `hf` client owns Hub access and downloads. BioV verifies the complete
selected bytes, adds portable companions and publishes a new directory without
replacing an existing one. `biov model inspect` verifies the saved bundle offline.

This is a bounded Rust CLI/library slice. It does not run models, infer scientific
meaning, resolve mutable revisions or claim that the selection is a complete or
usable model. Initial download publication supports Linux x86_64. Inspect has no
download-platform gate, but other platform distributions and clients require
their own acceptance. The implementation and validation scope are SPEC D16/T74.

## Download explicit files

Choose the full 40-hex Git commit and every exact repository-relative filename.
Options precede the repository; there is no implicit whole-repository selection,
branch/tag resolution, shortened commit, wildcard expansion or recursive folder
download. The destination must not exist, unless it already contains the exact
verified BioV resource requested. Even a pre-existing empty directory is not
adopted. Quote filenames and directory names in your shell when necessary.

The small public fixture below selects configuration only:

```sh
biov model download \
  --revision f171d7baecaf37b5da5a3616d8833b9969753535 \
  --local-dir ./tiny-bert-config \
  hf-internal-testing/tiny-random-bert \
  config.json tokenizer_config.json

biov model inspect ./tiny-bert-config
```

This exact two-file download was exercised through the actual BioV CLI with
the official pinned fallback on Linux x86_64; see the acceptance evidence below.
Neither weights nor a complete tokenizer, license or model
card is selected. A successful configuration-file download establishes none of
those files' presence and does not establish inference readiness. Add any wanted
native documentation or other files explicitly in a **new** selection/directory.
Never infer model architecture, training data, biological versions, licensing or
fitness from a repository name or suffix.

Success writes one JSON object to stdout:

- `status`: `downloaded`, `reused` or `verified`
- `local_dir`: the current absolute resource directory on this execution host
- `resource`: the complete format-1 portable record described below

The host path is a convenience, not a durable identity or a remote file transfer.
Progress, upstream output, help and diagnostics use stderr. Errors return status
2 with no success JSON; diagnostics retain an upstream failure's status rather
than inventing success. This route does not expose a model MCP tool or a Python
model-loading API.

## Use the upstream client

BioV prefers an installed `hf` on `PATH`. Compatibility requires a stable
three-integer version at least 0.34.0 and below 3.0.0, and `hf download --help`
must expose `--revision`, `--repo-type` and `--local-dir`. BioV queries the actual
client version and records it. An explicitly selected missing or incompatible
client is an error, never permission to substitute another executable.

```sh
biov model download \
  --revision f171d7baecaf37b5da5a3616d8833b9969753535 \
  --local-dir ./tiny-bert-installed-hf \
  --hf /path/to/hf --no-install \
  hf-internal-testing/tiny-random-bert \
  config.json tokenizer_config.json
```

Replace `/path/to/hf` with an existing executable. `--hf FILE` takes precedence
over `BIOV_HF_BIN`; absent both, `PATH` discovery applies. Relative executable
paths are resolved in the caller's original working directory. The child retains
that working directory, so relative upstream cache/token/tool paths retain their
usual meanings.

If no compatible implicit `hf` is available, BioV uses an already installed uv:

```text
uv tool run --no-config --no-python-downloads \
  --from huggingface-hub==2.1.1 --with 'httpx2[socks]' hf ...
```

uv owns this on-demand environment and its cache. BioV requires the fallback to
report version 2.1.1; it does not copy a client, create an installation registry,
download Python, modify the user's installed `hf` or edit shell settings. uv and
a suitable existing Python must be available. The primary package is pinned;
its remaining dependencies use uv's upstream resolution, not a BioV full
transitive lock. The added upstream `httpx2[socks]` extra supplies SOCKS transport
support; it does not change the pinned primary client. `--uv FILE` overrides
`BIOV_UV_BIN`, otherwise `uv` is found on
`PATH`. BioV does not bootstrap a missing uv for this route.

`--no-install` disables fallback provisioning. It is **not** an offline-download
flag: a new resource still invokes the installed `hf`, with that client's
ordinary network/authentication/cache behavior. Existing exact resource reuse
and `inspect` need no `hf`, uv, Python or network access at all.

The actual acquisition command targets a fresh sibling staging directory:

```text
hf download REPO FILE... --repo-type model --revision COMMIT \
  --local-dir /absolute/fresh/staging-directory
```

BioV adds no token argument and stores no credential in its record. Existing
Hugging Face authentication, environment variables, access restrictions,
transport and native cache metadata remain the upstream client's responsibility.
Configure account access separately with official Hugging Face tools when needed.
There is no BioV login, gated-access approval, upload or custom Hub endpoint
interface. A custom `HF_ENDPOINT` or an enabled `HUGGINGFACE_CO_STAGING`
(`1`, `ON`, `YES` or `TRUE`, case-insensitive) is rejected for new downloads;
the initial source contract is `https://huggingface.co`. Official-client cache reuse does not
become an independent upstream authenticity check.

For upstream usage, see the [official Hugging Face CLI guide](https://huggingface.co/docs/huggingface_hub/guides/cli)
and [uv tool guide](https://docs.astral.sh/uv/guides/tools/). BioV's deliberately
narrow command grammar above is not a pass-through for every `hf` option.

## Portable files and provenance

```text
tiny-bert-config/
  config.json                         # unchanged selected native bytes
  tokenizer_config.json               # unchanged selected native bytes
  BIOV_MODEL_RESOURCE.json             # relative inventory and acquisition facts
  BIOV_MODEL_RESOURCE_README.md        # scope, limits and standalone reader
  .cache/huggingface/                  # optional upstream private metadata
```

Provider-relative paths are preserved, including nested paths. The record and
README are companion files; they do not replace provider metadata. Preserve both
companions and every selected payload together when copying or moving a bundle.
The optional `.cache/huggingface` subtree is excluded from the verified inventory,
is not traversed by verification and is unnecessary for reading or exact BioV
reuse. Other unrecorded payload files are errors. Symlinks, special files and
missing/changed selected files are rejected.

`BIOV_MODEL_RESOURCE.json` format version 1 records:

- `provider: "hugging_face"`, `repository_type: "model"`, repository and exact
  caller-supplied revision
- `revision_provenance: "caller_supplied_full_git_commit_passed_to_hf"` and
  `scope: "selected_files"`
- Sorted unique `selected_files`, with one `inventory` entry per file containing
  relative `path`, complete `bytes` and lowercase `sha256`
- `acquisition`: actual `hf` client/version, successful exact-file/exact-revision
  fresh-directory invocation observed by BioV, no BioV payload transformation,
  local companion creation time and unknown (`null`) upstream download time
- `metadata`: model meaning not interpreted, selected native metadata
  authoritative when present and no model code execution
- `verification`: SHA-256 coverage of every selected file's complete bytes,
  authenticity not established, scientific QC not performed, no independent
  revision resolution and the excluded upstream metadata subtree
- `companions`: companion filenames, README byte count and README SHA-256

The schema is strict: unknown fields, changed scope/acquisition declarations and
unsupported versions are not accepted. The format-1 README is an immutable
companion template, not a scratch document to edit in place. Its checksum is
validated, and native `inspect` requires the supported template's exact bytes.
Notes added as extra payload files also invalidate the exact inventory; keep
personal annotations outside the saved bundle.

The creation timestamp means **when the local companion was recorded**. It does
not become a provider publication time or the original download time. Software,
record-schema and biological versions are separate facts. The supplied commit
was passed to `hf`; BioV does not independently resolve or authenticate the
provider's returned revision. Checksums establish consistency with the supplied
record, not producer authenticity, a signed supply chain, licensing, scientific
quality or model safety. A changed record can describe changed bytes.

## Move and read without BioV

Copy or move the full resource directory to another location. The companion
README contains a complete standard-library Python reader that validates the
bounded format-1 schema, README identity, exact relative inventory and every
selected file's complete bytes. It then summarizes all files by native filename
suffix and total bytes. Copy its displayed shell block, or save its displayed
Python as `verify_model_resource.py` outside the bundle and run:

```sh
python3 -I -S verify_model_resource.py ./moved-tiny-bert-config
```

This reader needs neither BioV nor Hugging Face nor the original download path,
catalog or network access. It streams files without decoding or loading weights.
An opaque weight format remains opaque; checksum verification does not license
unpickling, importing downloaded Python or enabling remote model code.

For the two-file example, after running the complete companion verifier, the
following ordinary JSON reader reports the native configuration fields actually
present. It reads both complete JSON documents and filters the configuration to
a few directly named metadata fields without loading a model:

```sh
python3 -I -S - ./moved-tiny-bert-config <<'PY'
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
names = ('config.json', 'tokenizer_config.json')
documents = {}
for name in names:
    with (root / name).open('rb') as stream:
        raw = stream.read(1024 * 1024 + 1)
    if len(raw) > 1024 * 1024:
        raise ValueError('this small-fixture JSON example has a 1 MiB file bound')
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ValueError('expected a native JSON object: ' + name)
    documents[name] = value
config = documents['config.json']
print(json.dumps({
    'native_document_key_counts': {name: len(value)
                                   for name, value in documents.items()},
    'selected_native_config_fields': {key: config[key]
        for key in ('model_type', 'hidden_size', 'num_hidden_layers', 'vocab_size')
        if key in config},
}, sort_keys=True))
PY
```

These values are retained native metadata, not independently verified claims
about trained weights, biological meaning or execution dependencies. The JSON
example is deliberately fixture-specific; it is not a reader for arbitrary
selected model formats. Its 1 MiB decode bound is separate from the generic
streamed payload verifier. Run it only on the already verified unchanged bundle
in a trusted directory.

## Publication, reuse and failure

New downloads use a fresh `.biov-model-download-*` directory beside the final
destination. Destination ancestors must be direct directories, without symlink
aliases. Only a successful upstream invocation with every exact selected
regular file present can produce companions and a published resource. Publication
uses Linux atomic no-replace rename. If another writer creates the destination,
that directory is preserved; this attempt fails and retains its unpublished
staging directory rather than adopting the winner.

An existing successful resource is reusable only after offline complete
verification and exact equality of repository, revision string and sorted file
selection. Reuse neither rewrites the record nor contacts upstream. A different
selection, damaged bundle or unrelated directory fails unchanged. There is no
refresh, repair, overwrite, model update, remove or cache-cleanup command.

Failure after staging starts retains the unpublished files and reports their
location on stderr. When no successful resource record was written, an
ordinary `BIOV_INCOMPLETE_MODEL_DOWNLOAD.json` diagnostic may describe the
incomplete attempt. This is not a successful resource or a portable acquisition
record. Review or retry native files with official `hf` as needed; a later BioV
attempt uses a fresh directory. BioV does not silently fall back to stale local
payloads after an upstream failure. Failed directories are not automatically
garbage-collected or published, and no SIGKILL recovery guarantee is implied.

## Bounds and acceptance scope

- 1–1,024 explicit selected files; each relative path is at most 512 UTF-8 bytes,
  with at most 32 components and at most 256 KiB aggregate selected path bytes
- Repository identifiers are bounded ASCII names or namespace/name identifiers,
  at most 192 bytes total and at most 96 bytes per component
- No absolute/traversing/unnormalized paths, glob syntax, duplicate paths,
  file/directory collisions, selected `.cache` files or reserved companion names
- JSON record at most 1 MiB, README at most 32 KiB, at most 65,536 traversed
  directory entries outside excluded upstream metadata
- Client version/help query stdout at most 128 KiB per query
- Complete-byte hashing uses a 128 KiB buffer. No selected payload-size or
  disk-consumption cap is implemented, and metadata bounds are not a measured
  whole-process memory or speed guarantee
- Filesystem checks require a trusted root without hostile concurrent writers;
  they are not an OS sandbox or content scanner

Source tests cover the exact selection/record contract, moved standalone-reader
verification, corruption/path/schema rejection and thin CLI delegation with a
fake upstream client. Test presence is not proof that actual `hf` networking or
installation worked. On 2026-10-05, the fresh Linux full-workspace Rust test run
passed, including 10 model unit tests and 12 model CLI cases. The actual BioV
CLI's uv fallback with official `huggingface-hub==2.1.1` and `httpx2[socks]`
downloaded the exact public fixture above successfully. Its complete selected
payload identities were:

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `config.json` | 548 | `54f9f0001ec62b6a53072a7c9d1b630e346ae51fca14a5799e54f0d918281872` |
| `tokenizer_config.json` | 321 | `9a8ed9b01c8a56b555dcdc31cd526ba9e488cfc16ee5a101b3ce64af81d34f3d` |

These 869 selected bytes identify that explicit configuration-only acquisition.
Both `tests/test_model_resources.py` cases passed with the source-built native
Linux binary and the real-public-download opt-in enabled. The real bundle was
moved, verified and reused offline under libseccomp network-call denial; isolated
`python -I -S` with BioV unavailable verified both complete payloads, summarized
all 869 bytes and read the complete native configuration/tokenizer JSON. The
source-built workspace also passed Clippy with warnings denied. Strict MkDocs
and all nine pre-commit hooks passed. The Linux wheel/sdist boundary checks and
release-wheel rebuild from the frozen source distribution passed; the unpacked
sdist's complete Rust suite passed (241 tests, including the doctest; three
explicit opt-in real-tool tests ignored). Final source/sdist Rust runs used one
test thread after intermittent `ETXTBSY` in unchanged temporary-executable router
fixtures on this cloud filesystem; an earlier ordinary parallel full run passed.

Both the independently installed production wheel and the rebuilt-sdist wheel
passed 507 installed acceptance cases outside the checkout, with one separate
real-Pixi lifecycle case explicitly skipped. These runs include both model
resource cases with actual public acquisition, native sequence/metrics and paired
Python/CLI behavior. The complete Python suite passed on 3.12, 3.13 and 3.14:
918 passed and 28 explicit opt-in skips on each. The latter two runs used the
installed ABI3 wheel. The separate real-Pixi lifecycle gate, required CI and
heterogeneous external review are not claimed as completed locally. Other
platform releases remain unvalidated. No model inference test or general
large-model/network/platform compatibility claim follows from this tiny fixture.

`tests/test_model_resources.py` is the separate installed-CLI/moved-reader
acceptance gate. `BIOV_TEST_BINARY` must name the native executable under test;
without it these checks skip explicitly. The default case uses a fake client;
`BIOV_TEST_REAL_MODELS=1` separately enables the small public download. Optional
`BIOV_TEST_REAL_UV` selects the official uv executable. For example, from a
checkout with its test dependencies installed:

```sh
BIOV_TEST_BINARY=/path/to/installed/biov \
  uv run --locked pytest -q tests/test_model_resources.py

BIOV_TEST_BINARY=/path/to/installed/biov BIOV_TEST_REAL_MODELS=1 \
  uv run --locked pytest -q tests/test_model_resources.py
```

The test moves the bundle so its original location disappears, denies network
calls in the offline child, clears manager lookup and runs the companion reader
with isolated standard-library Python and BioV absent. The opt-in real case also
reads both complete native JSON documents without importing a model library.
Passing a source-built binary validates that binary; it does not by itself
validate an independently installed wheel or sdist.
