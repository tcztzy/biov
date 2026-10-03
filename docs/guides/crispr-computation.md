# Deterministic CRISPR computation

BioV now distributes the former GEEPilot `crisprprimer` package, its original
scoring and restriction-enzyme lookup tables, the existing Docker report bridge,
and the Azimuth scoring runner.
GEEPilot owns the task skills: safety screening, selection of biological questions
and methods, interpretation, and follow-up analysis. This is a software ownership
migration, not a new validation of the scientific models or presets.

## Interfaces

- `crisprprimer` and `python -m crisprprimer`: existing Python workflow and CLI.
- `crisprprimer.nuclease`, `crisprprimer.score`, and
  `crisprprimer.id_converter`: existing deterministic Python interfaces.
- `crisprprimer-docker` and `python -m crisprprimer.docker`: existing Docker
  invocation and HTML-to-JSON/CSV conversion; native arguments remain available
  through `--help`.
- `biov-azimuth` and `python -m biov.azimuth`: existing 30-base input validation,
  Docker/fixed-fork selection, batch scoring and JSON/CSV/TSV output.

The library requires Python 3.12+. The Azimuth runner invokes the preexisting
legacy image or its explicitly pinned Python 3.10 fork environment; those
scientific environments are separate from BioV's Python runtime. The Docker
recipe is packaged as `biov/assets/azimuth.Dockerfile` and can be located with
`importlib.resources.files("biov").joinpath("assets/azimuth.Dockerfile")`.
Existing backend availability and model-pickle compatibility limits still apply.

The unpublished NAU-to-MSU mapping is intentionally not bundled: its source and
redistribution provenance have not been verified. `NAU2MSU` and `MSU2NAU` globals
are not exported. RAP/MSU conversion continues to use the published RAP-DB
mapping; NAU identifiers are not supported. The design workflow accepts MSU IDs
or genomic regions; RAP annotation lookup was not implemented in the migrated
workflow and now reports that limitation explicitly.

## Execution and cache ownership

The Python workflow invokes native BLAT with `biov.software.run_software` and
reads its native headerless PSL file. It retains its previous BLAT settings and
propagates a nonzero exit before parsing or caching output. BioV's configured
Pixi manifest must support BLAT on the execution platform; the bundled production
manifest currently targets Linux. The adapter is covered with native-output
fixtures; this migration does not claim a new live BLAT validation on macOS.

This API owns local temporary input/output files and rejects configured SSH
native-command dispatch. For another execution host, run the complete script on
that host using the existing BioV execution interfaces. No local temporary path
is implicitly copied to another host.

BioV and native fsspec settings own download caches. Importing crisprprimer no
longer overwrites the application's fsspec cache. Azimuth uses `BIOV_HOME/azimuth`
by default or the explicit `--cache-root` argument; the former
`GEEPILOT_AZIMUTH_CACHE_ROOT` variable is no longer used. Caller-specified BLAT
cache directories retain their existing format and must be kept separate for
different reference inputs.

## Verification

`tests/test_crisprprimer.py` checks public imports, bundled resources, CLI metadata,
original fixed score values, the native BLAT arguments and PSL parsing, compressed
reference input, cache reuse, process failure propagation and the local-file
boundary. `tests/test_azimuth.py` retains the existing input/backend-selection
regressions. Packaging CI checks both Python packages and their assets.

These checks do not execute a new design experiment or demonstrate measured guide
efficiency or specificity. The existing two-step managed sequence analysis remains
covered separately by `tests/test_analysis_example.py` and its opt-in real MCP
acceptance tests.

## Paired-read alignment

`biov.align_paired_reads(reference_fasta, fastq_r1, fastq_r2, bam_path)` owns the
shared BWA-to-BAM execution used by GEEPilot's Cas9202602 analyses. It builds the
reference index and runs BWA-MEM through `biov.software.run_software`, then uses
pysam to sort and index the BAM. `mem_args` and `sort_args` preserve the caller's
native parameters; use `build_index=False` with query-name sorting. Native BWA
stderr remains visible, and command or pysam failures propagate.

BWA's `-o` option writes the temporary SAM directly; its behavior is defined in
[the locked BWA 0.7.19 source](https://github.com/lh3/bwa/blob/v0.7.19/fastmap.c#L180).
The default coordinate-sorted result includes a BAM index. The API requires local
paths and rejects SSH command dispatch; the packaged BWA environment currently
supports Linux only. Running the complete script on a suitable host is separate
from sending local filenames to a remote command.

`tests/test_alignment.py` supplies explicit SAM fixtures at the native execution
boundary and checks actual pysam sorting/indexing, arguments and failure handling.
These checks do not claim a live BWA run or a reanalysis of experimental data.
