---
name: single-cell-annotation
description: Analyze single-cell RNA-seq data using measured counts, integration, reference transfer, and cell-type annotation. Use for scRNA-seq QC, embeddings, marker interpretation, reference mapping, or annotation confidence.
license: CC BY 4.0
---

# Single-Cell Cell Type Annotation

Three complementary approaches; combine them rather than trusting one.

1. **Manual marker-based.** Score clusters against known marker genes
   (Scanpy/Seurat). Useful for small datasets and novel cell types; confidence
   depends on marker specificity and supporting evidence.
2. **Automated classifiers.** Pre-trained tools assign labels directly. Fast
   and consistent, but bounded by their training data.
3. **Reference-based mapping.** Map query cells onto an annotated reference
   atlas. Best when a good reference exists for the tissue.

Practical rules:

- Annotate at cluster level first; refine ambiguous clusters at cell level.
- Never rely on a single marker gene — require a panel and check specificity
  against all other clusters.
- Report confidence and flag clusters that remain ambiguous instead of forcing
  labels.

Read `references/single_cell_annotation.md` for the full best-practices guide
(distilled from Luecken & Theis et al., *sc-best-practices.org*), including
tool choices, marker databases, and evaluation strategies.

## Counts, integration and transfer

Keep cells and genes uniquely identified in AnnData. Retain original counts in
a named layer, record organism, assay, gene namespace, donor, batch and QC
exclusions, and avoid densifying large sparse matrices. Assess doublets,
ambient RNA and low-quality cells using the actual experiment. Integration
must retain biological differences; batch and condition confounding cannot be
resolved by a prettier UMAP.

| Task | Direct method and inputs | Output / checks |
| --- | --- | --- |
| scVI / scANVI | Use scvi-tools with raw counts, declared batch and, for scANVI, label/unlabelled-category keys. Fit on the intended training/reference data. | H5AD with latent coordinates and cell IDs; retain training/convergence metadata and uncertainty. Do not train count likelihoods on log-normalized values or invent labels for unlabelled cells. |
| Harmony | Use harmonypy on the stated PCA representation with explicit batch covariates. | Corrected embedding with unchanged cell order; assess batch mixing and biological conservation separately. Harmony coordinates do not replace measured expression for differential testing. |
| Reference transfer / PopV | Supply annotated reference, query, common gene namespace, batch keys and selected methods. Use PopV's native preprocessing/training and retain per-method predictions. | Labels, agreement/confidence, unknowns and reference version. Query cells must not leak into held-out reference validation. |
| Pan-Human annotation | Use the official panhumanpy model/reference and feature-name requirements. | Hierarchical labels, confidence and reference coverage. Verify organism/tissue compatibility. |
| UCE / IMA mapping | Use the official UCE preprocessing, vocabulary and pinned checkpoint; compare to a reference generated with the same model/preprocessing. | H5AD embedding and neighbour distances/votes with cell IDs. Equal embedding dimensions alone do not make different models comparable. Out-of-reference cells can remain unknown. |
| STATE | Use the official State/SE checkpoint and documented count/gene inputs. | Embedding with cell IDs, checkpoint, model vocabulary and preprocessing record. Do not repeatedly halve batches on failure or silently alter gene IDs. |
| TranscriptFormer | Supply the matching species/model, gene IDs, raw counts and official vocabulary handling; explicitly choose cell/gene embedding and layer. | H5AD/numeric embedding, excluded/duplicate gene report and model revision. Do not synthesize an assay label, assume transformed counts are raw or drop duplicate genes silently. |

Save H5AD plus a JSON record of selected parameters, input references and
artifact paths. Foundation-model embedding generation is prediction, not
evidence that its cell labels are correct; inspect marker consistency and
reference coverage before interpretation. Models and reference atlases are
external assets and need explicit versions, resource budgets and license checks.

Native APIs: [Scanpy](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html),
[scvi-tools](https://scvi-tools.org/),
[PopV](https://github.com/YosefLab/popV),
[panhumanpy](https://github.com/satijalab/panhumanpy),
[UCE](https://github.com/snap-stanford/UCE),
[State](https://github.com/ArcInstitute/state), and
[TranscriptFormer](https://github.com/czi-ai/transcriptformer).

Migration source for these workflow decisions: Biomni commit
`400c1f366b96a35ca253e13c9b06c5076af41d65`, `biomni/tool/genomics.py`.
These instructions replace orchestration; no old function signature or
model performance claim is preserved.

Run Scanpy, scVI/scANVI and Harmony scripts with
`biov exec scvi -- annotation.py`. Use `biov exec popv -- transfer.py` or
`biov exec panhumanpy -- annotation.py` for those annotation methods.
For UCE, STATE or TranscriptFormer, invoke the corresponding `uce`, `state` or
`transcriptformer` tool through `biov exec TOOL -- ...` with its native arguments.
UCE needs existing absolute paths for `--adata_path`, `--model_loc`,
`--spec_chrom_csv_path`, `--token_file`, `--protein_embeddings_dir` and
`--offset_pkl_path`; specify `--dir` as an absolute output directory ending in `/`.
Use the files supplied with the selected model revision.

Read and write files in the working directory or configured data directories.
For a specialized tool without BioV, first obtain input files with BioV in a
`biov run prepare.py` script in the BioV application environment, then pass their paths to the tool.
