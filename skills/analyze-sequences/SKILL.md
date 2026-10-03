---
name: analyze-sequences
description: Analyze DNA, RNA, protein, and plasmid sequences with established libraries and command-line tools. Use for sequence annotation, primer and assembly checks, editing outcomes, motifs, structural features, conservation, or phylogeny.
---

# Sequence analysis

Read the actual sequences and their assembly/accession versions first. Keep
strand, circularity, genetic code, and coordinate convention explicit. Use
BioV for file access and interval operations; call the analysis library itself.
Save sequences as FASTA/GenBank, annotations as TSV/Parquet, structures as
mmCIF/PDB, and a small JSON summary with methods, versions and artifact paths.
Do not serialize Python objects or describe an in silico result as an experiment.

## Annotation, alignment, and editing

| Task | Inputs and execution | Results and checks |
| --- | --- | --- |
| Retrieve CDS or plasmid sequence | Resolve species and accession/version before fetching with NCBI Entrez or the depositor's sequence download. A gene symbol can have multiple isoforms; an Addgene numeric ID is not an NCBI accession. Parse CDS features with Biopython rather than assuming the first search hit. | FASTA plus accession, transcript/protein IDs and selected CDS feature. Preserve alternative CDSs unless the user chooses one; verify translation with the annotated genetic code. |
| ORFs | Use Biopython sequence translation in the requested frames and genetic code; state start-codon and minimum-length rules. Report reverse-strand intervals in the original sequence coordinate system. | ORF coordinates, frame, strand, nucleotide and peptide sequence. Overlapping/nested ORFs are separate candidates, not confirmed genes. |
| Plasmid / bacterial annotation | Use pLannotate for plasmids and Prokka for bacterial assemblies. Supply circularity, organism metadata and database versions; run the native CLI with checked exit status. | GenBank/GFF3, nucleotide/protein FASTA and annotation tables. Keep the tool's confidence and database provenance. |
| Short-sequence matching / variants | Use `Bio.Align.PairwiseAligner` for scored alignment and an appropriate read aligner for sequencing data. Set mismatch/gap scores explicitly and search both strands when appropriate. | Aligned sequences and substitutions/insertions/deletions in reference coordinates. A position-by-position `zip` cannot detect indels; report ambiguous equal-score placements. |
| CRISPR outcomes | Supply amplicon reference, nuclease, guide/PAM orientation and sequenced reads; use CRISPResso2 for read-level quantification. For two already assembled sequences, inspect their alignment and the actual cut interval. | Allele counts, denominators, indel/substitution coordinates and read QC. A donor's short matching substring does not establish HDR, and a cut-site match does not predict delivery or editing efficiency. Use the existing `sgrna-design` skill for guide selection. |
| Comparative genomes / haplotypes | Align assemblies or reads to a versioned reference with the appropriate aligner; use a real variant caller and phasing method when haplotypes are requested. | Alignment, VCF and phase evidence. Do not compare unaligned sequence indexes or infer ancestry/divergence time from variant counts alone. |
| Barcode reads | Parse FASTQ with Biopython; require the barcode regex or both flanks, orientation, expected length and quality threshold. Count exact observed barcodes before applying an explicitly chosen error-correction method. | Read and barcode counts, excluded reads, correction mapping and abundance table. Edit-distance clusters alone are not cell lineages. |

## Primer and assembly checks

- **Primer design and verification coverage:** use Primer3 with the supplied
  template, allowed/excluded regions, product size, salt and primer concentrations.
  Retain Tm, GC, hairpin and dimer results. Map both primers to the full template
  and check off-target products against the relevant reference. For Sanger
  verification, draw the directional coverage intervals using the user or
  instrument's usable read length; do not assume 800 usable bases for every read.
- **PCR simulation:** use pydna with both primer sequences and explicit circularity.
  Report every valid product and its sequence; missing primer binding is a failed
  simulation, not permission to return the requested target region. A simulated
  band image is not gel evidence.
- **Restriction maps and digestion:** use `Bio.Restriction.RestrictionBatch`
  and enzyme metadata with explicit linear/circular topology. Preserve cut
  offsets, overhang sequences and fragment lengths; verify that fragment lengths
  sum to the template length.
- **Golden Gate oligos and assembly:** use the enzyme's actual Type IIS cut
  offsets and DnaCauldron's `Type2sRestrictionAssembly`. Enumerate backbone
  and insert cuts, internal sites, overhang compatibility and orientation. Check
  the assembled sequence, junctions, expected copy counts and circular closure;
  ambiguous compatible overhangs are multiple possible products. Protocol
  selection belongs in `lab-protocols`, not a hard-coded reaction recipe.
- **Codon redesign:** use DNA Chisel with the correct genetic code, host codon
  usage table/version and requested sequence constraints. Verify unchanged
  translation, frame and stop codon. Maximizing each codon's frequency is not
  evidence of increased expression.
- **Provided genome edits:** represent user-specified sequence edits against one
  reference coordinate system and apply them without mutating the input plan.
  Save the edited FASTA/GenBank and coordinate mapping. The old therapeutic-delivery
  function changed positions twice and claimed transformation readiness; neither
  that algorithm nor the biological efficacy claim is retained.

## RNA and protein

| Task | Method | What to retain / avoid |
| --- | --- | --- |
| RNA folding and structure features | ViennaRNA `RNA.fold_compound`, MFE and partition-function methods; parse pair tables for paired/unpaired counts and stems. | Sequence, temperature, parameter set, dot-bracket, energy and ensemble uncertainty. Match sequence/structure lengths; do not mix unsupported pseudoknot notation with a nested-only parser or sum arbitrary base-pair energies. |
| Protein conservation and phylogeny | MAFFT for MSA; IQ-TREE for model-based phylogeny; Biopython for alignment/tree I/O. | Alignment, gap-aware conservation definition, tree/model and support values. Never replace a failed MSA with padded strings or a collection of unrelated pairwise alignments. |
| Structure comparison | Match chain sequences/residues, then use `Bio.PDB.Superimposer` on the corresponding atoms. | Atom correspondence, coverage, transformed structure, RMSD and per-residue displacement. Matching residue numbers alone is insufficient across unrelated numbering systems. |
| Disorder | Use IUPred's documented model and mode or its official service. | Per-residue scores, threshold and contiguous intervals. Preserve sequence indexing; disorder predictions are not measured structure. Check the model's license before use. |
| N-glycosylation sequons | Scan all overlapping `N-X-S/T` windows with `X != P`. | Positions and matched residues, explicitly labelled sequons rather than occupied glycosites. |
| O-glycosylation | Use a stated, versioned predictor such as NetOGlyc with its applicable organism and license; a local S/T fraction may be reported only as sequence composition. | Per-residue score and model threshold. The upstream arbitrary window/proline heuristic is not a glycosylation probability. |
| FAS / other protein domains | Run HMMER against a versioned profile database, with protein sequence and significance thresholds. Translate nucleotide input only with a justified frame/code. | Domain accession, residue interval, score/E-value and overlap handling. A substring match to a domain name does not establish enzymatic activity. |
| ChatNT questions | Invoke the pinned `InstaDeepAI/ChatNT` model using its documented Transformers interface when this model is specifically needed. Inspect any required remote code before executing it. | Model revision, input sequence and generated answer; label it as a model response, not a validated sequence annotation. |
| Protein embeddings | Use the official ESM model/featurizer with an explicit layer, checkpoint, sequence-length policy and pooling excluding special/padding tokens. Preserve isoform IDs before any documented gene-level aggregation. | Numeric NPZ/HDF5 arrays plus accession/isoform index and model metadata. Report excluded/truncated sequences; do not silently skip failed proteins or use pickle-based tensor dumps as the analysis interchange format. |

Use the official [Biopython documentation](https://biopython.org/docs/latest/Tutorial/index.html),
[Primer3](https://primer3.org/), [pydna](https://github.com/pydna-group/pydna),
[DnaCauldron](https://edinburgh-genome-foundry.github.io/DnaCauldron/ref/user_classes.html),
[ViennaRNA](https://www.tbi.univie.ac.at/RNA/), and
[DNA Chisel](https://edinburgh-genome-foundry.github.io/DnaChisel/) for executable APIs.
Specialist predictions require the selected model and reference data; they
are not function-compatible replacements for the old tools.

Migration source: Biomni commit `400c1f366b96a35ca253e13c9b06c5076af41d65`,
`molecular_biology.py`, `genetics.py`, `biochemistry.py`, `glycoengineering.py`,
and the sequence-related functions in `genomics.py`, `microbiology.py`,
`synthetic_biology.py`, `biophysics.py`, and `systems_biology.py`.

Run sequence analysis scripts with `biov exec python -- analysis.py`.
Use `biov exec esm -- embedding.py` or `biov exec chatnt -- query.py` for
those selected models. Invoke native sequence tools with `biov exec TOOL -- ...`
and their documented arguments.
Available native tools are `plannotate`, `prokka`, `crispresso2`, `mafft`,
`iqtree` and `bwa`. The `hmmer` entry already runs `hmmscan`, so pass only its
arguments: `biov exec hmmer -- profiles.hmm proteins.fasta`. Select another
HMMER program with native Pixi arguments and the same manifest, for example
`biov pixi run --manifest-path "$BIOV_ENVIRONMENT_MANIFEST" --environment hmmer -- hmmsearch proteins.hmm proteins.fasta`.

Read and write files in the working directory or configured data directories.
For a specialized tool without BioV, first obtain input files with BioV in a
`biov run prepare.py` script in the BioV application environment, then pass their paths to the tool.
