# biov-storage

Bounded offline registration and resolution of complete **native file trees**.
This crate has no Polars, Tokio, network client, database or downloader. It copies
an existing directory into a portable immutable snapshot; it does not convert
formats or infer missing scientific metadata. CLI/MCP adapters use these public
Rust operations:

```rust,ignore
let store = biov_storage::NativeStore::new(store_root)?;
let registered = store.register(source_root, request)?;
let resolution = store.resolve(resolve_request)?;
```

Both roots must already exist. Construction and resolution are read-only.
Registration alone creates store layout. `source_path` is a normal relative path
beneath the configured source root, which must be disjoint from the store root.
Original files are not changed, moved, deleted, or hardlinked. Copied permissions,
ownership, timestamps, extended attributes and sparse allocation are not part of
the native-byte contract. Ordinary readers may use complete copied files directly.

## Portable layout

```text
artifacts/<namespace>/<canonical-accession>/snapshots/sha256-<64hex>/
  README.md           # direct analysis entry-file list and ordinary-reader example
  acquisition.json    # operational receipt, schema version 1
  checksums.sha256    # ordinary sha256sum -c checksums.sha256
  source/             # complete original tree, including native README/empty dirs
.staging/             # unfinished work is never ready
.locks/               # cooperating writer locks
```

Native README/catalog/records remain authoritative. Wrapper files do not replace
or rewrite them. A complete snapshot remains understandable and readable when
moved without BioV, source directories, a catalog database or network access.
There are no derived indices in this slice. Receipt paths are source-relative;
response paths are store-relative plus a clearly labeled execution-host path.
Original source paths are not stored or required for resolution. The receipt
includes `readme_sha256` for the wrapper entry point, so verification does not
depend on recreating its Markdown with a future BioV version.

`Registration` and `Resolution::Ready` contain compact `SnapshotSummary` values,
not entire inventories. `Resolution` also represents `Miss`, `Ambiguous`,
`Unavailable`, and `Corrupt`; errors distinguish invalid requests/packages,
resource limits, filesystem failures, declaration conflicts and unsupported
publication. There is no silent latest choice. Unversioned GCF references may
match multiple biological versions, and one biological version may have multiple
annotation/package snapshots. A caller may pin the full source snapshot ID.

`requested_ref` in a receipt is the first registration request. The registration
response separately preserves the current request when verified content is reused.
`registered_at_unix_seconds` is local registration time. Unknown acquisition
facts (`acquired_at`, `source_url`, `acquisition_tool`) are null. Local copying
does not establish how, when or from whom the input was downloaded. Neither
checksums nor a receipt authenticate the provider or prove scientific correctness.
Resolution repeats native checks and verifies the stable `validation.method`
identifier under the receipt schema. Historical `validation.limits` prose is
bounded, untrusted description, not an integrity/authentication claim; editorial
changes to that text do not invalidate unchanged native content. A changed method
is rejected rather than permitting claims of a different validation contract.

## Exact source identity encoding

Inventory entries include all files and directories below `source/`, but not
`source/` itself. Paths use UTF-8 and `/`, sorted by UTF-8 bytes. Empty directories
participate in identity. Wrapper metadata, registration time, filesystem metadata
and original root paths do not participate. The SHA-256 input is:

1. ASCII `biov-native-tree-v1` followed by one NUL byte
2. For each sorted entry: one ASCII `D` or `F`; unsigned 32-bit big-endian path
   byte length; raw UTF-8 path bytes; unsigned 64-bit big-endian file byte count;
   32 raw SHA-256 bytes
3. Directory byte counts are zero and their 32 digest bytes are all zero

Receipt directory entries have `kind: "directory"`, `bytes: 0`, `sha256: null`;
file entries have `kind: "file"` and lowercase SHA-256 hex. Hashes cover saved
bytes, including compressed `.gz` bytes, not decoded biological content.

A source digest identifies source content, not declared scientific semantics.
Identical bytes under one canonical accession with a changed scope or changed
representation mapping, **including reordered paths**, produce an explicit
`DeclarationConflict`. Existing receipts are never silently overwritten or
reinterpreted. Mappings preserve native catalog/caller vector order. Different
source content creates a different snapshot; no in-place annotation refresh exists.

## Native adapters and limits

RefSeq accepts an extracted NCBI Datasets V2 package containing native README,
`md5sum.txt` and `ncbi_dataset/data/dataset_catalog.json`. Every accession-bearing
catalog group must match the explicit versioned canonical GCF reference. Report
only groups containing `DATA_REPORT` may omit accession; an accession-bearing
group must contain files. Catalog paths are relative to
`ncbi_dataset/data`, not package root. Every catalog file must exist, be nonempty,
have the exact declared byte length, and be covered by supplied MD5 checksums;
the catalog itself must also be MD5-covered. All supplied MD5 targets are checked.
Dehydrated or incomplete native packages are rejected. Native file types map to
`genome_fasta`, `rna_fasta`, `cds_fasta`, `protein_fasta`, `gff3`, `assembly_report`
and `sequence_report`; other safe types become `native_<lowercase-type>`.
Missing RNA means representation unavailable, not absence of RNA genes. No
sequence or annotation scientific QC is performed.

PDB accepts classic four-character entry accession syntax, canonicalized to
uppercase, plus caller-declared `entry` or `assembly:<positive decimal ID>` scope
and a mapping from representation names to ordered vectors of source-relative
files. Files must exist and be nonempty. This deliberately does **not** parse or
verify mmCIF/PDB identity, assembly coordinates, format, chain semantics, provider
origin, or compressed/decoded equivalence. Extended PDB accessions, entity/chain
identity and other providers are outside this slice. Resolve defaults to RefSeq
`assembly` and PDB `entry`; an assembly must be explicitly selected.

Native source bytes have no universal size cap. Copying and hashes use 64 KiB
buffers; registration performs a second full source inventory to detect ordinary
in-flight modifications/additions/removals, and resolution fully rehashes source
files and repeats native validation. Reuse is also fully verified. These are
currently linear full-byte operations, not a large-data performance claim.

Bounds are 100,000 source entries; 4,096 bytes per relative path; 16 MiB per
native metadata/receipt/wrapper file; 128 representation names; 256 paths per
representation; 1,024 namespace entries for unversioned discovery; 1,024 snapshot
directory entries per unpinned discovery; 64 ambiguity candidates; and 64 KiB
serialized response. Retained inventory metadata is additionally charged at path
bytes plus 256 bytes per entry, with a 16 MiB budget before copying more files.
Versioned GCF/PDB and explicit snapshot pins use direct
paths, so unrelated entries do not impede exact lookups. Exceeding bounds returns
a structured error, never truncated scientific data. Metadata is still retained
in memory and these bounds are not a universal peak-RAM guarantee.

## Filesystem and publication assumptions

All paths must be normal, bounded UTF-8 relative paths. Traversal, absolute
paths, control characters, backslashes, Windows-reserved path syntax/device
names, symlinks and special files are rejected. Both configured roots and every
selected path are checked. Complete trees include all regular files, even those
not catalogued; native metadata remains untouched. Trusted root permissions are
required: these checks are not a sandbox against hostile concurrent writers.

A unique staging directory is completed and validated first. Publication uses an
established `fs4` exclusive per-canonical-accession lock and an atomic no-replace
directory rename (Rustix `NOREPLACE` on Linux/Android/Apple, non-replacing directory
rename on Windows). Existing snapshots are verified and either reused or rejected;
corrupt winners are never replaced. Immutability is an API publication contract,
not operating-system write protection; later edits to ordinary snapshot files are
detected by verification. Other targets return unsupported publication.
Files, directory trees, newly created parent entries and publication parents are
synced on Unix. Linux behavior is tested; other platforms/filesystems, network
filesystems, power-loss recovery and hostile-writer races are not validated.
Process-crash remnants in staging are ignored, not automatically deleted. There
is no claim of a cross-filesystem transaction or comprehensive crash recovery.
