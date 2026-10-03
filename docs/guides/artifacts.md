# Identifier-backed analysis

BioV's computation boundary is deliberately small: it maps a persistent
identifier to a validated file in the environment running the code. It does not
wrap GC, alignment, parsing, or other analysis functions. The LLM can therefore
write ordinary Biopython or command-line code it already knows, while BioV hides
cache, download, and executor-local path details.

The packaged `artifact_capabilities.json` is the source of supported namespaces,
artifact kinds, defaults, and file formats. Call `biov.artifact_capabilities()`
to inspect it without network or cache access. It describes fixed file
representations, not arbitrary API requests.

These formats minimize avoidable download and conversion work. They do not
make every dataset immediately suitable for every analysis: XML and SOFT still
need their domain parsers, and experiment-specific normalization and biological
interpretation belong to the workflow. BioV does not invent missing PDFs,
matrices, structures, or experimental file selections.

## Read through fsspec

Installing BioV registers the manifest namespaces as read-only fsspec
protocols. Use the complete identifier URI; importing BioV first is not
required. Every protocol uses the same provider files and cache as `biov.path`
and `biov.open`.

```python
import fsspec

with fsspec.open("uniprot://P42212", "rt") as sequence_file:
    fasta_text = sequence_file.read()

with fsspec.open("uniprot://P42212", "rt", artifact="entry_json") as entry_file:
    entry_json = entry_file.read()
```

Pass `artifact` as a filesystem option to select a non-default file. The
existing BioV readers also accept these URIs:

```python
import biov

records = biov.read_fasta("uniprot://P42212")
annotation = biov.read_gff3(
    "refseq.gcf://GCF_000006945.2",
    storage_options={"artifact": "annotation_gff3"},
)
```

`read_fasta` keeps its existing dictionary of Biopython sequence records;
`read_gff3` returns a `BioDataFrame`. File objects support `read` and `seek`.
`fsspec.open_local(uri)` returns the real cached filename, and
`fs.info(uri)` returns the URI, byte size, and `type="file"`. Metadata lookup
can fetch the selected artifact on a cache miss. Write modes are rejected
before downloading; directory listing is unsupported. Provider and identifier
errors propagate unchanged.

For example, `biov.path("pubmed://23193287")` returns a PDF filename and
`biov.open("dbsnp://rs328", mode="rt")` opens full RefSNP JSON. MCP discovery
for these namespaces returns a JSON file description without downloading it;
the same URI passed to BioV or fsspec reads the actual data.

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

The same cache can also hold complete entry metadata, fetched on its own first request:

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
`Artifact.path` points to the original member that `dataset_catalog.json`
selects for the requested kind, such as the `GENOMIC_NUCLEOTIDE_FASTA` file for
`genome_fasta` or the `GFF3` file for `annotation_gff3`;
`Artifact.package_root` points to the package root. It also exposes the
requested and catalog-canonical identifiers, kind, and file size. One cached
package serves every RefSeq kind without another download.

A valid package-directory cache hit reads the local catalog and file metadata,
then skips `datasets`. A cache miss downloads and extracts in a sibling staging
directory and publishes the whole package atomically. When
`BIOV_MAX_FILE_BYTES` is set, the archive and every extracted member must stay
within it; an oversized package is rejected before anything is published. An
explicit version such as `GCF_000006945.2` must match the catalog exactly. For a
versionless request, the single versioned assembly returned in the official
package becomes the canonical identifier.

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
`pdb://1EMA`. Request `artifact="alphafold_cif"` on a UniProt URI for a predicted
model; zero or multiple AlphaFold models raise an error. BioV does not choose
among them or substitute an experimental structure.

## NCBI and other individual files

PubMed resolves the PMID with the official
[PMC ID Converter](https://pmc.ncbi.nlm.nih.gov/tools/id-converter-api/), selects
the current version, and reads its
[PMC cloud metadata](https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/). Only its
published PDF object is downloaded, with PMID, version, and MD5 verification.
The former OA web service was discontinued in August 2026 and is not used.
A PMID without an available PDF fails; no publisher page is scraped as a
fallback.

ClinVar uses [EFetch VCV XML](https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/)
for a numeric Variation ID. dbSNP uses the complete
[RefSNP JSON service](https://api.ncbi.nlm.nih.gov/variation/v0/). Neither
returns the old ESummary subset.

GEO preserves the decompressed native
[Series Matrix or SOFT file](https://www.ncbi.nlm.nih.gov/geo/info/download.html).
The default requires exactly one nonempty matrix for a GSE accession. Multiple
matrices require selecting an explicit upstream file; GSM/GPL/GDS accessions
require `artifact="soft"`. BioV does not infer a study's required platform,
join matrices, or select supplementary count files.

Reactome retrieves SBML from
`https://reactome.org/ContentService/exporter/event/{accession}.sbml`.

Individual files are staged, validated, and published atomically under
`$BIOV_HOME/artifacts/<namespace>/<accession>/`. A valid cache hit skips the
provider request. HTTP and network failures propagate. A cached file is a
snapshot; an unversioned accession can change upstream. Record the file,
accession, retrieval date, and any relevant provider version with the analysis.

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

For a local locked Pixi environment, MCP `run_analysis` can instead execute an
ordinary Python script and preserve its inputs, outputs, checks and execution
record. See [managed analysis](analysis.md) for that request contract and bounded
previews. Identifier discovery remains read-only; analysis execution is explicit.
Managed analysis does not dispatch through SSH or LSF.
