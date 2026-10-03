---
name: evaluate-drug-candidates
description: Evaluate compounds using cheminformatics, docking, prediction models, measured pharmacology data, and sourced drug safety records. Use for descriptors, poses, ADMET, binding predictions, drug repurposing, stability, release, or preclinical response analysis.
---

# Drug candidate research

Record compound identity/stereochemistry, target identity, model or assay
version, units and the task's comparison. Separate measured endpoints, model
predictions and evidence lookups. Use native scientific tools, keep their
structured/tabular outputs, and return JSON with provenance and artifact paths.
This skill does not supply patient-specific prescribing or certify safety.

## Structures and predictions

| Task | Inputs / direct tool | Output and checks |
| --- | --- | --- |
| Molecular properties | Parse and sanitize SMILES with RDKit; explicitly choose fragment/salt, charge, tautomer and stereochemistry handling before calculating descriptors. | Canonical identity and named descriptors/units (MW, cLogP, TPSA, HBD/HBA, rotatable bonds). Do not silently neutralize a molecule; Lipinski descriptors are not measured ADMET. |
| Binding-site discovery | Use the selected AutoSite release on a prepared receptor; inspect pockets and the coordinate system. | Pocket coordinates/volumes and score with units. A proposed box is a candidate site, not established ligand binding. |
| AutoDock Vina | Prepare receptor and ligand with Meeko, retaining protonation, cofactors and stereochemistry; supply box center/size in Å, exhaustiveness, seed and CPU count. | PDBQT/SDF poses, docking scores and run metadata. Check a relevant reference ligand when available; docking score is not measured affinity. |
| DiffDock | Use the official versioned model/checkpoint and valid receptor plus ligand representation; retain preprocessing and sampling parameters. | Poses and model confidence with structural QC. Do not claim confidence ranks are binding energies. |
| ADMET / sequence-based affinity | Use the chosen DeepPurpose model with its matching featurizer, endpoint, training data and checkpoint. Confirm SMILES/protein validity and applicability. | Per-compound predictions with endpoint units and model version. Do not round, rescale or call a score a probability without that model's contract. |
| TxGNN repurposing | Use the official model or a verified exported prediction table with disease/drug IDs and model/data version. | Ranked scores and provenance; exact disease mapping must be inspectable. Do not treat a sigmoid-transformed score as calibrated efficacy. |

Follow the official [Vina preparation/docking workflow](https://autodock-vina.readthedocs.io/en/stable/docking_basic.html),
[RDKit](https://www.rdkit.org/docs/GettingStartedInPython.html),
[DiffDock](https://github.com/gcorso/DiffDock),
[DeepPurpose](https://github.com/kexinhuang12345/DeepPurpose), and
[TxGNN](https://github.com/mims-harvard/TxGNN).
Verify weight/data terms separately from code licenses and select only
the model required by the task.

## Measured pharmacology

- **Drug release:** supply measured concentrations, sampled/replaced volumes,
  initial loading and times. Convert to cumulative mass with sampling correction
  before calculating fraction released. Fit justified release models on their
  applicable domains with uncertainty/residuals. The observed maximum is not
  automatically 100% loading; an early-time fit must not be interpreted across
  the full release curve.
- **Accelerated stability:** fit measured potency/degradation and physical
  endpoints at the actual temperatures/humidity. Extrapolation needs a justified
  kinetic/Arrhenius model and validation. The upstream hard-coded stability
  curves are retired; formulation names alone cannot provide shelf life.
- **Xenograft growth:** use longitudinal per-animal measurements with named
  treatment/control and baseline. Choose a repeated-measures/mixed model or
  prespecified endpoint; report the exact TGI formula, attrition and uncertainty.
  The first alphabetic group is not an implicit control; repeated observations
  are not independent animals.
- **Radiolabeled biodistribution:** retain time, tissue mass, administered
  activity, decay-correction convention and biological replicates. Report
  %ID/g or the actual measured unit, tumor-to-organ ratios and a justified
  time–activity fit/AUC. Infer clearance only if volume/activity scaling permits it.
- **Dosimetry:** require calibrated time–activity curves, physical decay,
  integration/extrapolation convention and radionuclide/geometry-specific
  S-values. Compute time-integrated activity and absorbed dose with explicit
  units; do not substitute arbitrary S-factors or turn this research calculation
  into a clinical dose recommendation.
- **VCOG grading:** obtain the applicable published VCOG-CTCAE version and
  endpoint-specific thresholds, units, reference intervals and clinical context.
  Save the source criterion with each assigned grade; ungradable observations
  stay ungraded. The upstream incomplete numeric dictionary and generic
  mild/moderate mapping are not the standard.
- **Chondrogenic aggregate assay planning:** use a cited method matching the
  cell source and assay. Report the supplied experimental plan or analyze its
  measured endpoints; do not claim to have cultured cells or generated results.
  Use `lab-protocols` for available reference procedures.

For blots use `analyze-biomedical-images`; for CpG associations use `analyze-omics`.

## Drug evidence and interaction records

Use exact database drug IDs and documented synonym mapping, retaining unmatched
and ambiguous names. Read a verified DDInter export as tabular/JSON data, preserve
the actual interaction text, severity, mechanism, management and source version.
No record means unknown, not “safe”. Candidate alternatives need the clinical
indication and evidence; membership in a broad class and absence of an entry do
not make a drug a safe replacement.

For FDA events, labels and recalls, use the
[official openFDA API](https://open.fda.gov/apis/) with explicit endpoint,
query/date fields, filters and pagination, and preserve source JSON. Apply date
filters to the actual query (the old function only printed some requested dates).
Record totals and truncation. Report adverse-event counts as reports, not
incidence or comparative risk. A disproportionality calculation needs a defined
full comparison table, duplicate handling and uncertainty; an arbitrary threshold
on a few drug-specific reports is not a validated safety signal.

Migration source: Biomni commit `400c1f366b96a35ca253e13c9b06c5076af41d65`,
`pharmacology.py` and the release-analysis function in `bioengineering.py`.
These are task instructions and direct ecosystem use, not BioV compatibility APIs.

Run RDKit and measured-data statistics scripts with
`biov exec python -- analysis.py`. Invoke the selected docking or prediction
tool with `biov exec TOOL -- ...` and its documented arguments.
The `vina-meeko` entry task invokes Vina, so use
`biov exec vina-meeko -- ...` for docking. For Meeko, set the same project manifest
used during setup and select its native program with
`biov pixi run --manifest-path "$BIOV_ENVIRONMENT_MANIFEST" --environment vina-meeko -- mk_prepare_ligand.py ...`.
Use `autosite` for pocket detection and
`biov exec deeppurpose predictions.py` for a complete prediction script. Library
entry tasks run Python; inspect the imported package rather than interpreting
Python help as package validation. The [shipped manifest](https://github.com/tcztzy/biov/blob/main/src/biov/assets/environments/pyproject.toml)
declares entry commands. DiffDock uses
`biov exec diffdock -- --config /absolute/path/to/inference.yaml ...`;
the YAML must name the selected model directories with absolute paths and
contains the authoritative inference settings.

Read and write files in the working directory or configured data directories.
For a specialized tool without BioV, first obtain input files with BioV in a
`biov run prepare.py` script in the BioV application environment, then pass their paths to the tool.
