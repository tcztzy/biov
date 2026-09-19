# Identifier-backed analysis

BioV's computation boundary is deliberately small: it maps a persistent
identifier to a validated file in the environment running the code. It does not
wrap GC, alignment, parsing, or other analysis functions. The LLM can therefore
write ordinary Biopython or command-line code it already knows, while BioV hides
cache, download, and executor-local path details.

The provider pairs are:

```text
refseq.gcf × genome_fasta → NCBI Datasets genomic FASTA
uniprot   × protein_fasta → UniProtKB entry FASTA
uniprot   × entry_json    → complete UniProtKB entry JSON
```

Other data types should be added as namespace × artifact providers, not as one
BioV function per biological computation.

## Use normal Biopython code

`path` accepts the same exact forms recognized by the MCP parser. GCF
forms include `refseq.gcf://GCF_000006945.2`,
`refseq.gcf:GCF_000006945.2`, and the allowlisted bare
`GCF_000006945.2`. UniProt accessions remain explicit, for example
`uniprot://P42212` or `uniprot:P42212`; bare `P42212` is not inferred. The
returned `Artifact` implements `os.PathLike`, so libraries see an ordinary file:

```python
import biov
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

genome = biov.path("GCF_000006945.2", artifact="genome_fasta")

gc_weighted_bases = 0.0
total_bases = 0
for record in SeqIO.parse(genome, "fasta"):
    gc_weighted_bases += gc_fraction(record.seq) * len(record.seq)
    total_bases += len(record.seq)

print(gc_weighted_bases / total_bases)
```

The artifact kind defaults from the namespace. Ordinary protein code is just
as direct:

```python
import biov
from Bio import SeqIO

gfp = SeqIO.read(biov.path("uniprot://P42212"), "fasta")
print(gfp.id, len(gfp.seq))
```

The same cache exposes complete entry metadata without another request:

```python
import json

import biov

with biov.open("uniprot://P42212", artifact="entry_json", mode="rt") as entry_file:
    entry = json.load(entry_file)

pdb_ids = [
    reference["id"]
    for reference in entry["uniProtKBCrossReferences"]
    if reference["database"] == "PDB"
]
```

Only the added `biov.path` call is BioV-specific. The model remains free to
use Biopython, pysam, samtools, or other established tools for the actual
analysis. `biov.open(..., mode="rt")` is available when a library wants a
file handle instead of a path.

`parse_identifier` is the strict, network-free single-value parser. Unlike the
MCP `parse_identifiers` prompt tool, it rejects surrounding prose and multiple
IDs.

## Official NCBI package and cache behavior

The RefSeq provider invokes the official NCBI Datasets CLI directly:

```console
datasets download genome accession GCF_000006945.2 \
  --include gff3,rna,cds,protein,genome,seq-report
```

BioV adds only `--filename` to place the ZIP in a temporary cache directory and
`--no-progressbar` to keep program output clean. It does not replace this with a
REST downloader.

After validating ZIP paths, BioV extracts the complete package without
flattening, renaming, copying, or adding a manifest. For example:

```text
$BIOV_HOME/artifacts/refseq.gcf/GCF_000006945.2/
├── README.md
├── md5sum.txt
└── ncbi_dataset/
    └── data/
        ├── assembly_data_report.jsonl
        ├── dataset_catalog.json
        └── GCF_000006945.2/
            ├── *_genomic.fna
            ├── genomic.gff
            ├── rna.fna
            ├── cds_from_genomic.fna
            ├── protein.faa
            └── sequence_report.jsonl
```

The files actually supplied by NCBI depend on the assembly. The returned
`Artifact.path` points to the original `GENOMIC_NUCLEOTIDE_FASTA` member named
by `dataset_catalog.json`; `Artifact.package_root` points to the package root.
It also exposes the requested and catalog-canonical identifiers, kind, and file
size.

A valid package-directory cache hit reads the local catalog and file metadata,
then skips `datasets`. A cache miss downloads and extracts in a sibling staging
directory and publishes the whole package atomically. An explicit version such
as `GCF_000006945.2` must match the catalog exactly. For a versionless request,
the single versioned assembly returned in the official package becomes the
canonical identifier.

Set `BIOV_HOME` to choose the base directory; otherwise BioV uses its normal
platform cache directory. The [`datasets` CLI must be installed](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
and available on `PATH` in the environment that runs the script. A missing or
rejected CLI command and a malformed package use separate stable artifact
errors.

## Official UniProt entry and cache behavior

The UniProt provider uses the
[documented individual-entry REST URL](https://www.uniprot.org/help/api_retrieve_entries)
for both original representations:

```text
GET https://rest.uniprot.org/uniprotkb/P42212.fasta
GET https://rest.uniprot.org/uniprotkb/P42212.json
```

Each representation is requested and cached only when needed. BioV validates
the exact requested accession before atomically publishing that individual
file:

```text
$BIOV_HOME/artifacts/uniprot/P42212/
├── P42212.fasta
└── P42212.json
```

The response byte streams and accession filenames are unchanged; BioV adds no
manifest, reformatted JSON, or derived cross-reference list. A JSON-only or
FASTA-only directory is a valid partial cache. Reading `entry_json` never
downloads FASTA, and reading `protein_fasta` never downloads JSON. Each existing
file is validated independently and skips its corresponding REST request.

A UniProt entry can cross-reference many experimental PDB structures with
different mutations, ligands, chains, and conditions. BioV therefore does not
silently choose one or place an invented `P42212.pdb` beside the sequence.
Experimental coordinates belong to an explicit PDB identifier such as
`pdb://1EMA`; a one-to-one predicted model belongs to an AlphaFold DB artifact.

## Run the complete script in its target environment

Run locally:

```console
biov run analysis.py -- GCF_000006945.2
```

The script inherits the current environment, working directory, and standard
streams. Its exit code becomes the `biov run` exit code.

Submit the same script to LSF:

```console
BIOV_LSF_PYTHON=/shared/venvs/biov/bin/python \
  biov run --executor lsf --queue short --job-name gc-analysis \
  --stdout 'logs/gc.%J.out' --stderr 'logs/gc.%J.err' \
  analysis.py -- GCF_000006945.2
```

BioV passes the interpreter, absolute script path, and arguments directly to
`bsub`; it does not build a shell command. It pins the submission working
directory with `-cwd`. A successful command returns a numeric LSF job ID and
means only that submission was accepted. The job may still be pending, running,
or later fail; inspect it with the site's normal LSF tooling.

The script path, working directory, selected Python environment, installed
libraries, required downloader/network access, and any configured `$BIOV_HOME`
must be visible where the job runs.
Artifact path resolution happens on the execution host. Consequently, a compute
node path is not returned unless the analysis code emits it. Sites with pre-populated
shared storage can add a provider for the same namespace × artifact contract
without changing generated Biopython analysis.

## LLM workflow

The intended division of work is:

1. MCP `parse_identifiers` discovers and validates IDs in the user's prompt.
2. The LLM writes ordinary analysis code, retaining the ID and adding
   `biov.path` at the file boundary.
3. `biov run` executes the complete script in local or LSF context.

The MCP server does not execute generated code. This separation keeps resource
discovery read-only while deterministic execution remains explicit and
auditable.
