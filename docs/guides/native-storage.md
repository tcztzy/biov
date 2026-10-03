# Offline native snapshot storage

BioV's storage rule is simple: copy the saved bundle to another machine and use
ordinary scientific tools to understand and analyze it without BioV. Provider
files, metadata and package layouts remain authoritative. The filesystem is the
durable record; BioV adds discovery, integrity checks and only the operational
facts missing from the native package.

## Implemented slice and boundaries

The `biov-storage` Rust library is a transport-independent, local copy-registration
and offline-resolution slice. It accepts already-present native source directories
under a configured source root and publishes immutable snapshots under a separate
configured store root. It does not download data, unpack remote archives, run an
HTTP service or replace the existing Python artifact cache. Thin `biov-rs`
CLI/MCP adapters expose the same Rust storage contracts. Calling the library is explicit; it does not discover and import existing
caches automatically.

The semantic adapter is deliberately narrow:

- **RefSeq:** an exact versioned `refseq.gcf` accession, resolved through an
  existing NCBI Datasets package's native catalog, with native checksums and local
  content checks
- **PDB declaration-only registration:** caller-declared identity and
  representation-to-file mappings for a structurally different native package. This
  preserves files and verifies their local integrity. It is not a PDB semantic
  validator and does not establish the claimed entry, chain or assembly from
  coordinates
- **UniProt, AlphaFold and GEO:** the real inspections below inform the design;
  their semantic bundle adapters and cross-record validation remain future work

This is server-side filesystem plumbing, not evidence of a deployed distributed
storage system. No durable database or catalog is introduced. Discovery scans
snapshot records on demand. Moving a store and resolving from its new root tests
reconstructed filesystem discovery; there is no database to delete or rebuild.

Source files are copied, never moved, deleted or hardlinked. Registration must
not rewrite custom local metadata or native biological records. Existing Python
cache roots, configured environment roots and analysis-output roots remain
separate. There is no automatic `imports/` migration or external-path registrar.

## Identity, layout and independent use

A readable namespace/accession directory groups immutable snapshots. The complete
native package lives beneath `source/`, where its own README, checksums, metadata
and filenames cannot collide with wrapper files. The wrapper README lists analysis entry files by representation, points to native
metadata and explains the declared scope, so a reader can start without crawling
the catalog. It does not create a converted analytical view. Snapshot references
use relative file paths; original absolute locations are not required to read a moved bundle.

```text
<configured-store-root>/
  artifacts/
    refseq.gcf/GCF_000005845.2/snapshots/sha256-<64-hex-digest>/
      acquisition.json
      checksums.sha256
      README.md
      source/
        README.md
        md5sum.txt
        ncbi_dataset/data/
          dataset_catalog.json
          assembly_data_report.jsonl
          GCF_000005845.2/
            GCF_000005845.2_ASM584v2_genomic.fna
            genomic.gff
            protein.faa
            cds_from_genomic.fna
            sequence_report.jsonl
```

This is an illustrative native payload from the inspected RefSeq acquisition,
not a promise that every package has those representations. Directory and receipt
schema details are defined by the library contract below. Native catalog paths
are relative to `ncbi_dataset/data`, not the source root. A catalog's report-only
group need not have an assembly accession.

Keep these identities distinct:

1. The biological namespace and canonical accession, including its own version
2. The requested representation and scientific scope
3. The exact native-source inventory identity and individual file hashes
4. Registration time and any separately known upstream acquisition time
5. Native provider release, annotation, entry, sequence or model version evidence

An assembly version does not freeze all annotation bytes. A PDB entry ID does not
freeze its later revisions. A local registration timestamp is not an upstream
release or retrieval date. Missing acquisition facts stay unknown; copying an
old package today must not imply it was downloaded today.

File preservation is distinct from format decoding or preview. This slice's
PDB declarations validate local bytes and paths, not whether an arbitrary file
can safely be parsed. A `ready` resolution means the declared files are locally
materialized and verified against the receipt; it does not establish scientific
reader support or biological validity. Future generic imports may retain opaque
files, including serialized formats such as pickle, with explicit unsupported
decoding status. Untrusted pickles must never be automatically loaded or executed.
Prepared Arrow tables are optional, not a universal replacement for FASTA,
BAM/CRAM, structures or other native formats.

The receipt is an operational inventory, not a universal biological object. It
points to native metadata and records identity, local file paths, sizes, hashes,
registration and validation scope. It must not duplicate sequences, annotations,
features or matrices, and must distinguish caller declarations from checks BioV
actually performed. Checksums establish consistency, not provider authenticity,
trusted provenance or scientific quality.

## Registration, offline resolution and failure behavior

Registration copies into a private staging area on the publication filesystem,
validates the complete staged package and publishes only a ready snapshot. A
partial stage is never a cache hit. Identical concurrent registrations reuse a
verified published winner; different native bytes create distinct snapshots and
cannot overwrite an earlier one. Immutability is the library's write contract,
not a filesystem permission sandbox against an administrator modifying files.

An offline lookup resolves a requested representation to an ordinary local path.
No DNS, HTTP request or download is attempted. It never selects whichever snapshot
was most recently modified. Multiple matching snapshots need an explicit selector;
a missing representation is different from a missing identity or damaged local
data. No error may turn an older valid snapshot into a partial replacement.

Native RefSeq file availability follows the catalog and actual files. Requesting
RNA during an earlier download does not make `rna_fasta` present. A missing RNA
entry is not evidence that the genome has no RNA genes. Metadata-only/dehydrated
packages cannot be reported as ready sequence data.

The configured roots are trusted local directories. Path validation rejects
escapes; it does not protect against hostile concurrent local writers. Staging,
copying and hashing can consume disk even when an eventual registration fails.
No implicit eviction or garbage collection is provided. Ordinary readers can
inspect a published snapshot without keeping a BioV process alive.
The copy preserves file bytes, relative names and empty directories; it is not
a backup of filesystem ownership, permissions, timestamps or extended attributes.
Validation is currently Linux-specific. The path checks do not establish safe
movement onto every case-folding or Unicode-normalizing filesystem; such targets
need collision checks and their own platform acceptance.

## CLI and MCP adapters

Registration and resolution take a JSON request file matching the Rust request
fields below:

```sh
biov-rs storage register --store-root STORE --source-root SOURCES --request-file register.json
biov-rs storage resolve --store-root STORE --request-file resolve.json
```

`STORE` and `SOURCES` must already exist and be disjoint trusted directories.
Request `source_path` is relative to the configured source root. Resolving a
moved store needs only its new store root and a resolution request. Request JSON
is limited to 64 KiB. Successful results are compact JSON on stdout; malformed
requests and operational errors use stderr and nonzero exit status. Resolution
`miss`, `ambiguous`, `unavailable` and `corrupt` are normal structured results,
not CLI transport failures.

The existing native MCP startup accepts optional `--store-root DIR` alongside
`--data-root DIR --output-root DIR`. Its `storage_register` and `storage_resolve`
tools use the same request/response contracts. MCP registration confines sources
to `--data-root`; the tool input cannot change that root. Both storage tools are
advertised even without `--store-root` and then return an explicit not-configured
tool error. These are native Rust adapters;
the existing Python `biov mcp` and provider-download interfaces remain separate.
An execution-host path returned by resolution is not a transfer to a remote client.

## Rust API and on-disk contract

The library entry points are `NativeStore::new(store_root)`,
`store.register(source_root, RegisterRequest)` and `store.resolve(ResolveRequest)`.
The store root exists independently of the source root; reopening a moved store
needs only its new root. `source_root` is supplied for registration and must be
an existing trusted directory disjoint from the store. No new environment
variable or TOML setting is wired into the Python configuration in this slice.

Registration request fields:

- `source_path`: UTF-8 path relative to `source_root`, naming the complete package
- `requested_ref`: the user's explicit reference, retained separately from the
  canonical identity
- `canonical_ref`: the explicit exact canonical reference
- `declaration`: `{"provider":"refseq"}`, or `{"provider":"pdb","scope":"entry",
  "representations":{"structure_cif":["1crn.cif"]}}`. PDB also permits
  `assembly:<positive decimal ID>` scope. The paths are source-relative, and one
  representation may deliberately name multiple files

For RefSeq, use a canonical versioned reference such as
`refseq.gcf:GCF_000005845.2`. The native slice names the GFF3 representation
`gff3`; the legacy Python route uses `annotation_gff3`. Other recognized mappings
include `genome_fasta`, `rna_fasta`, `cds_fasta`, `protein_fasta`,
`assembly_report` and `sequence_report`; other safe catalog file types become
`native_<lowercase fileType>`. Native V2 catalog validation requires matching canonical GCF accession-bearing
groups, while permitting only `DATA_REPORT` members in accessionless groups.
An accession-bearing group must not be empty. Native README, MD5 file,
catalog, catalog member sizes and all declared catalog members must be present
and verified. This does not perform FASTA/GFF biological consistency checks or
infer an unavailable representation. PDB identity is restricted to classic four-character accessions, syntax-checked
and canonicalized to uppercase; declared files must exist, but claimed identity, coordinates, revisions, compression and
entry/assembly meaning are not parsed or independently checked.

Resolution fields are `reference`, `representation`, optional `snapshot_id` and
optional `scope`. Default scope is `assembly` for RefSeq and `entry` for PDB.
A ready response contains a compact snapshot summary and a vector of files:
`relative_path` is relative to the current store root; `execution_host_path` is
an ordinary absolute path on the resolving host. It is not a remote file transfer
or a durable reference to the original machine. Consumers may pass that path to
normal readers. The complete receipt stays on disk, rather than being returned
as an unbounded tool response.

Resolution reports explicit `ready`, `miss`, `ambiguous`, `unavailable` or
`corrupt` states. Ambiguity supplies bounded snapshot summaries; unavailable
includes the snapshot's available representations. Invalid input/package,
metadata/output limits, I/O, corruption, declaration conflict and unsupported
publication have distinct library errors. A supplied snapshot selector must be
an exact valid snapshot ID; no `latest` policy or mutable alias is implemented.

### Receipt and content identity

`acquisition.json` schema version 1 contains the canonical/requested references,
scope and declaration, native metadata entry points, representation mappings,
unavailable representations, inventory, validation method/limits and local
registration facts. The inventory has `path`, `kind`, `bytes` and `sha256` for
each source-relative entry. Directories have zero bytes and null hashes, including
empty directories. File hashes cover the exact copied bytes, including gzip when
retained. The original source location is not persisted. `acquired_at`,
`source_url` and `acquisition_tool` remain null because registration cannot verify
how an existing package was originally obtained.
`readme_sha256` separately identifies the exact wrapper README bytes; it is not
part of the native-source digest. Verification checks the stored entry point,
rather than regenerating it with a potentially newer documentation template.

`registered_at_unix_seconds` records the first successful local registration.
Reusing a snapshot does not rewrite that timestamp or its first `requested_ref`;
the registration response separately returns the current request and `reused`.
No new acquisition event or upstream check is implied by reuse.

The source identity hashes this deterministic binary encoding:

1. ASCII prefix `biov-native-tree-v1` followed by one NUL byte
2. All inventory entries sorted by UTF-8 path bytes, using slash-separated paths
3. For each entry: `D` or `F` byte; unsigned 32-bit big-endian path byte length;
   the UTF-8 path; unsigned 64-bit big-endian byte length; 32 raw SHA-256 bytes
4. Directory length and digest bytes are zero; the source root itself is excluded

The resulting lowercase SHA-256 is `source_content_sha256`; `snapshot_id` is
`sha256-` plus that full digest. Source directories and files are included, but
wrapper receipts, checksums and README are excluded. A short display prefix is
never the authoritative selector. Identical source trees with conflicting scope
or representation declarations fail with `DeclarationConflict`, preserving the
existing snapshot rather than silently changing its meaning.

`checksums.sha256` is a conventional file checksum inventory with paths under
`source/`. The receipt carries directory records and the complete content identity.
Source-native integrity records remain unchanged inside `source/`. This wrapper
is small relative to the payload and deliberately contains no biological record
reserialization.

### Limits, publication and verification cost

There is no aggregate scientific-file byte cap. Metadata and traversal are still
bounded: 100,000 source entries, 16 MiB per bounded metadata document, 4,096 bytes
per path, 128 representations, 256 resolved files, 1,024 scanned snapshot
candidates and 64 ambiguity candidates. Serialized public responses are limited
to 64 KiB. These are denial-of-service and output limits, not proof of peak-memory
usage or support for arbitrary provider package complexity.
Inventory traversal additionally charges path bytes plus 256 bytes per entry
against a 16 MiB metadata budget before retaining more entries. Exact canonical
references and snapshot pins use direct paths; bounded enumeration is required
only for unversioned-reference or snapshot discovery, and counts ignored entries.

Registration streams copies and hashes; resolution fully rehashes the snapshot
for verification. A lookup can therefore be expensive for a large native tree.
This is an explicit integrity trade-off for this bounded slice, not an optimized
large-store indexing/query design. Do not rehash a multi-terabyte input for each
future analytical batch; such work needs a separately validated verification and
index strategy.

The implementation uses a store-local filesystem lock and atomic no-replace
publication. A replacing rename is insufficient for immutable storage. Linux and
Apple publication use no-replace rename support, Windows uses its non-replacing
rename behavior, and unsupported platforms/filesystems return a dedicated error.
Only tested platforms are validated; these implementation branches are not
release or distributed-filesystem claims. Symlinks, special files, traversal and
control-character paths are rejected. Filesystem roots and concurrent writers
remain inside a trusted operator boundary.

### Analyze a moved RefSeq snapshot without BioV

Copy the complete snapshot directory, including the wrapper and native `source/`
tree, to an unrelated directory. From that copied snapshot:

```sh
sha256sum -c checksums.sha256
(cd source && md5sum -c md5sum.txt)
```

The native catalog can be read with Python's standard library to locate and
summarize every genome FASTA record. This example uses only the copied files:

```python
import json
from pathlib import Path

source = Path("source")
data = source / "ncbi_dataset" / "data"
catalog = json.loads((data / "dataset_catalog.json").read_text())
group = next(
    item for item in catalog["assemblies"] if item.get("accession") == "GCF_000005845.2"
)
member = next(
    item for item in group["files"] if item["fileType"] == "GENOMIC_NUCLEOTIDE_FASTA"
)
records = bases = gc = 0
with (data / member["filePath"]).open() as fasta:
    for line in fasta:
        if line.startswith(">"):
            records += 1
        else:
            sequence = line.strip().upper()
            bases += len(sequence)
            gc += sequence.count("G") + sequence.count("C")
print(
    {
        "records": records,
        "bases": bases,
        "literal_GC_fraction": gc / bases if bases else None,
    }
)
```

This simple summary counts literal G/C bases, not BioV's weighted ambiguous-base
GC definition. The larger existing
[RefSeq example](rust-migration.md#reproduce-the-real-refseq-package-inspection)
independently checks coordinates and the complete protein/CDS records. A copied
PDB mmCIF can be parsed directly with Biopython's `MMCIFParser`; that ordinary
reader's biological parsing is separate from declaration-only registration.

Acceptance must independently block network and remove BioV from the reader
environment, read complete records and perform a meaningful summary. Copying a
snapshot and checking only that filenames exist is not sufficient. See SPEC D13
for the bounded lifecycle scenarios and the validation status below.

## Verified implementation scope (2026-10-03)

The source-built Linux slice passes 32 storage-library tests and the full
130-test Rust workspace suite. Its 17 actual MCP/CLI subprocess cases pass
against the independently installed binary, including the dataset CSV correction.
The complete Python regression run with all local acceptance fixtures configured
passes 912 tests with two unrelated skips.

Four installed storage acceptance cases pass with kernel network syscalls blocked
for both the native CLI and independent readers: synthetic snapshot relocation and
revision selection; the actual RefSeq package; the actual PDB package; and a killed
registration of a 64 MiB + 1 byte file followed by fresh-process lookup and retry.
The tests preserve original fixture hashes, remove disposable import paths,
relocate the store, recompute the documented source inventory identity, check the
saved README hash, execute ordinary-reader examples, and read an individually
exported snapshot after its store is unavailable. No database is present.

The real RefSeq reader sees one 4,641,652-base genome, 4,300 proteins and 4,318 CDS
FASTA records. The PDB reader sees one model, chain A, 46 residues and 327 atoms.
These checks demonstrate the inspected packages, not generic scientific validation
by the registrar or coverage of all provider variants.

`tests/test_native_storage_portability.py` defaults to offline synthetic fixtures
when an installed `BIOV_TEST_BINARY` and independent `BIOV_STANDALONE_PYTHON` are
configured. Real cases additionally require `BIOV_TEST_REFSEQ_PACKAGE` and
`BIOV_TEST_PDB_PACKAGE`; PDB reading uses `BIOV_STANDALONE_BIOPYTHON` with Biopython
1.88 and no BioV. CI runs the synthetic cases without downloading provider data.
Explicit invalid configuration fails rather than silently skipping. The network
sandbox is confined to child processes and requires Linux/libseccomp; no host
security settings or privileges are changed. Other platforms, power-loss fault
injection and distributed filesystems remain outside the tested claim.

## Lessons from five real sources

These are measured observations of public files inspected on 2026-10-03, not
claims that all five providers are implemented or that these examples cover every
biological edge case. No downloaded biological payload is bundled with this guide.

### RefSeq: GCF_000005845.2

The official NCBI `datasets` CLI version 18.38.0 returned a 4,154,411-byte ZIP and
13,814,260 extracted bytes. Its native README, seven-entry `md5sum.txt`, catalog,
assembly report and sequence report already describe most of the package. Preserve
them rather than introduce a competing biological schema. The report uses
camelCase keys; CLI summary JSON is a different representation.

Independent full-file analysis found one 4,641,652-base genome, 4,300 protein
FASTA records and 4,318 CDS FASTA records. GFF sequence references and coordinates
matched the genome; the `thrL` interval `NC_000913.3:190..255` gave the expected
66-base CDS after explicit conversion from 1-based inclusive coordinates. GFF CDS
feature rows are not a count of distinct proteins. Requested RNA was absent from
the package and catalog.

The assembly release was 2013-09-26 and annotation release 2026-09-02. Therefore
`GCF_000005845.2` is not an exact annotation-content pin: changed annotation under
the same accession requires a different immutable snapshot. A moved package
produced identical independent analysis. That original inspection used no network
calls; it did not separately disable networking.

See [the reproducible RefSeq inspection](rust-migration.md#reproduce-the-real-refseq-package-inspection),
the [official assembly record](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/)
and [native package documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/).

### PDB: 1CRN

Native `1crn.cif.gz` and its losslessly decoded `1crn.cif` preserve PDBx/mmCIF
semantics. The observed entry revision was 1.5, dated 2024-10-30; Biopython read
one model, chain A, 46 residues and 327 atoms. Native metadata describes X-ray
diffraction and 1.5 Å resolution. The optional entry API JSON was pretty-printed
by the inspection client and must not be described as exact HTTP body bytes.
The mmCIF alone is sufficient to reconstruct the structure semantics.

Entry, polymer entity, chain instance and biological assembly are different
scopes. Preserve assembly IDs and distinguish mmCIF `label_asym_id` from author
chain identifiers. An entry coordinate file cannot masquerade as a requested
assembly file. RCSB, PDBe and PDBj are access points to the shared wwPDB archive,
not three independent experiments. This small example does not validate complex
assemblies, NMR ensembles, alternate conformations or RNA complexes.

The declared-package path preserves the caller's choices; biological verification
against mmCIF and compressed/decoded equivalence are not implemented semantic
checks. See the [official entry](https://www.rcsb.org/structure/1CRN),
[download and assembly services](https://www.rcsb.org/docs/programmatic-access/file-download-services)
and [PDBx/mmCIF dictionary](https://mmcif.wwpdb.org/).

### UniProt: P69905 — semantic adapter planned

Native FASTA and full entry JSON preserve distinct useful representations. The
observed entry had 142 amino acids, entry version 219, sequence version 2 and
release 2026_03. Both HBA1 and HBA2 occur in the record, so a gene-name key loses
identity. Preserve requested isoform accessions; sequence version 2 does not mean
isoform 2. Database release, entry version, sequence version and acquisition time
are separate fields.

A directory containing independently downloaded FASTA and JSON does not establish
that both represent a consistent provider observation. A future adapter must
cross-check their sequences/version evidence, preserve originals and report
uncertainty if upstream changed between requests. See the [native entry JSON](https://rest.uniprot.org/uniprotkb/P69905.json),
[FASTA](https://rest.uniprot.org/uniprotkb/P69905.fasta) and
[accession semantics](https://www.uniprot.org/help/accession_numbers).

### AlphaFold: AF-P69905-F1 — semantic adapter planned

The native prediction API array, v6 ModelCIF and v6 PAE JSON form a useful bundle.
The inspected model had 142 residues, 1,077 atoms and a 142 × 142 PAE matrix.
ModelCIF, API and UniProt sequences agreed for this observation. Resolve download
URLs from the API response rather than generate them from assumed filenames.

`F1` is a model fragment designation, not an isoform. Database model v6 is not the
reported `AlphaFold Monomer v2.0 pipeline` algorithm version. Explicit target
range and exact sequence are needed before joining a prediction to a current
protein record. The API returns an array even for one model; preserve one-to-many
results and do not replace requested isoforms with canonical proteins.

ModelCIF pLDDT is predicted local confidence, including when stored in atom B-factor
fields; it is not an experimental temperature factor or diffraction resolution.
PAE represents different uncertainty. The native `max_predicted_aligned_error`
field is not necessarily the observed maximum matrix value. An empty API array
is a negative observation, not a successful zero-byte model; transient errors
must not destroy an earlier usable model. These response/lifecycle semantics
remain future adapter work.

See the [prediction endpoint](https://alphafold.ebi.ac.uk/api/prediction/P69905),
[release notes](https://www.ebi.ac.uk/pdbe/news/alphafold-database-release-notes) and
[confidence interpretation](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/).

### GEO: GSE5623 — semantic adapter planned

The inspected native Series Matrix and family SOFT total 15,015,281 compressed
bytes. The matrix contains 22,810 probes × 24 GSM samples on GPL198. Family SOFT
already includes platform annotations and sample tables, so an extra GPL download
was unnecessary. All 547,440 matrix values matched the native sample values and
all probe IDs joined to platform IDs. Preserve native sample ordering, repeated
SOFT metadata keys, empty/control probes and multiple mappings: 979 annotation
rows contain multiple AGI mappings.

These are processed MAS5 microarray signals, not counts or raw reads. All 24
sample titles describe salt treatments; a broader study summary cannot establish
local control samples. Other series can have multiple platforms/matrices and
assay-specific supplements. Source record updates, platform annotation dates,
HTTP timestamps and acquisition time must remain distinct. No inference about
experimental design or normalization follows from successful storage.

See the [series record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5623),
[official files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5623/) and
[GEO download guidance](https://www.ncbi.nlm.nih.gov/geo/info/download.html).

### Integrity is representation-specific

The UniProt/AlphaFold inspection retained exact application-body bytes after
HTTP content decoding. `Content-Length` or provider storage hashes may describe
compressed wire/storage bytes instead; a mismatch with a decoded file is not
necessarily corruption. Original compressed downloads, decoded materializations
and reformatted JSON must not share an unqualified claim of byte identity.
A local copy-registration receipt cannot reconstruct unknown historical HTTP
headers, client options or download times.

## Analysis-ready derived layer: planned

Preserving native bytes is the starting point, not the full usability goal. The
current wrapper directly lists useful analysis entry files and their declared
scope. A future, separately identified derived layer should make complete tables,
indices and relationships easy to use while leaving native originals untouched.
No conversion, normalization or biological interpretation is implemented by the
current registrar.

Examples requiring their own contracts and independent tests:

- RefSeq: clearly related genome, annotation and reference sequence IDs, plus
  derived FASTA/tabix indices when appropriate and actually present
- GEO: complete expression, sample metadata and probe annotation views preserving
  GSM order, platform relationships, repeated metadata, MAS5 meaning and multiple
  probe-to-gene mappings; no fabricated counts or controls
- PDB: coordinates and native metadata with explicit entry versus assembly scope,
  assembly operations and author versus label chain IDs

The planned layers are acquired native data, prepared/materialized analysis-ready
views, and durable results. Each derived view needs exact native input hashes,
actual commands/parameters and code/tool versions, random seed when relevant,
output-affecting dependencies, ordered schema, row meaning, joins/shard order,
known units/coordinates/null semantics and an ordinary-reader example. A transform
key identifies the declared recipe and inputs; a separate output SHA-256 identifies
the actual result bytes. Uncaptured external state or nondeterminism prevents a
deterministic-reuse claim. Unknown facts stay unknown. Reuse raw snapshot references
rather than require a duplicate raw payload for every view. A portable export must
explicitly materialize the required dependency closure or declare an intentional
output-only scope; broken absolute paths or unavailable input hashes do not make
a self-contained reproducible analysis bundle.

### Hugging Face comparison: design evidence, not a dependency

Hugging Face distinguishes downloaded Hub objects from its processed Datasets
cache. Its Arrow-backed datasets demonstrate the utility of reusable prepared
views and memory-mapped access, but do not establish a BioV performance result or
justify converting every FASTA, BAM or structure into Arrow. See the official
[cache management](https://huggingface.co/docs/datasets/en/cache) and
[Arrow architecture](https://huggingface.co/docs/datasets/en/about_arrow) guides.

Datasets fingerprints track state and transformations rather than serving as a
cryptographic identity of complete output bytes. The inspected source at commit
`af4347bd438a00e7c688661f7ae1e377674e0a97` uses xxhash64 and incorporates state and
file modification times in relevant fingerprint paths. Keep recipe/cache keys
separate from source and output SHA-256 identities. This comparison was source
and documentation inspection; no installed Hugging Face experiment was run.
See [fingerprint source](https://github.com/huggingface/datasets/blob/af4347bd438a00e7c688661f7ae1e377674e0a97/src/datasets/fingerprint.py)
and the [cache concept guide](https://huggingface.co/docs/datasets/en/about_cache).

For BioV, future large-data interfaces should distinguish metadata inspection,
bounded scan, indexed access and explicit materialization under a memory budget.
Domain-native formats can retain their own useful access paths. Any exported
symlink/blob-based cache needs its referenced byte closure materialized before it
can count as portable. Future acceptance must move prepared bundles, read them
without BioV, and verify invalidation when input bytes, recipe or output-affecting
dependencies change, independently from output-integrity checks.

## Future slices, kept separate

- Transactional downloads through official provider tooling, safe archive handling,
  hydration, response evidence and explicit refresh; current registration never
  performs these operations
- Semantic UniProt/AlphaFold/GEO adapters, supported namespace decisions,
  paired-sequence/target-range checks and multi-platform/multi-model fixtures
- Managed local imports with explicit copied versus external-reference ownership,
  and result input-lineage integration; retain originals and mark unknown schema,
  units or biological meaning honestly
- Optional filesystem-rebuildable discovery indices and explicit selection/alias
  records; a future index must never be the sole identity or retention record
- Quotas, requested-download estimates and an auditable GC plan honoring owners,
  pins, active readers and retained input dependencies; no implicit data eviction
- Derived scientific indices, lazy/disk-backed execution and bounded batches,
  followed by profiling before remote/object-store or distributed features

Streaming a copy and hash is not large-table analytical support. The existing
Rust dataset route retains its independent 64 MiB retained-data charge and bounded
CSV/IPC parsing contracts. That charge is not a peak-memory guarantee. Storage
must not impose the same small universal limit on valid native genome files, but
it also must not promise constant whole-process memory, streaming scientific
parsing or distributed performance from a stream-copy implementation.

A future index should be rebuildable from durable files. Derived `.fai`, tabix,
BAM/CRAM or analytical indices need exact source digests and builder versions,
separate from immutable native payloads. Exact computational reproducibility
requires the input dependency closure or explicit unavailable-input declarations;
a digest or historical URL alone is insufficient.
