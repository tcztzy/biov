---
name: analyze-assay-data
description: Analyze measured biochemical, plate-reader, cytometry, growth, and chromatographic assays. Use for kinetic fitting, standards/calibration, flow gates, proliferation, titers, CFU, or quantitative assay comparisons.
---

# Analyze measured assays

Start from measured data, sample IDs, units, blanks/standards, biological
replicates and the requested comparison. Preserve raw data and exclusions.
Fit with SciPy/statsmodels or the assay's established package, save tabular
results and residuals, and return a JSON summary naming input/output files,
model, parameters, uncertainties and failed QC. Do not substitute simulations,
typical values or default calibration factors when observations are missing.

## Kinetics and calibration

| Assay | Inputs and method | Required checks / interpretation |
| --- | --- | --- |
| Protease / enzyme kinetics | Time-by-substrate measurements, fluorescence-to-product calibration and enzyme concentration in compatible units. Fit a justified initial linear window, then a bounded Michaelis–Menten or other stated model. | Vmax/Km with uncertainty and residuals. Compute kcat only from concentration/time rates and active enzyme concentration; raw RFU/s divided by enzyme concentration is not kcat. Modulator IC50 needs measured dose-response data. |
| ITC | Injection heats/volumes, cell volume, cell/syringe concentrations, temperature, dilution controls and instrument convention. Use a validated binding model with injection dilution/mass balance. | Kd, stoichiometry, ΔH, residuals and uncertainty; ΔG uses the thermodynamic standard-state convention and compatible energy units. Do not default unknown concentrations to 1 and 10 M or assume a 1.4 mL cell. |
| Circular dichroism | Wavelength/signal, baseline, concentration, path length and spectral units; optionally measured thermal curves. | Corrected spectra and fitted thermal parameters with reversibility/baseline checks. Protein secondary-structure fractions require a suitable reference-based deconvolution; a count of positive/negative wavelength bins is not a composition estimate. |
| ELISA / EBV antibody / ATP luminescence | Standards and samples with blanks, dilution factors and replicate IDs. Fit the validated linear or nonlinear response-to-concentration calibration over its working range. | Back-calculated standards, concentration and uncertainty; report values outside quantification limits. ATP amount additionally requires assay volume, then cell/protein normalization with explicit units. Antibody serostatus needs assay-specific validated cutoffs. |
| Crystal violet biofilm | Replicate absorbance, blank, negative/positive controls and sample design. Subtract blank and calculate the requested relative biomass. | Per-replicate and sample summaries. A one-sample t-test against a noisy control mean ignores control uncertainty; use the appropriate paired/group model. |
| OD growth curves | Time, blank-corrected OD, dilution and strain/replicate IDs. Select a measured exponential interval or fit an identified logistic/Gompertz model. | Growth rate, doubling time and model-specific lag with residuals/uncertainty. OD is not cell count without calibration; do not take logarithms of nonpositive corrected OD. |
| CFU dilution counts | Observed counts, actual dilution fraction, plated volume and replicate identities. Calculate count divided by dilution fraction and plated volume; combine countable plates with an explicit weighting rule. | CFU/mL, uncertainty, excluded/confluent plates and detection limit. Never generate Poisson counts from an assumed concentration and call them observed colonies. |
| HPLC–ICP-MS arsenic / GC fatty acids | Chromatograms or peak table, method-specific retention standards, response factors, internal standard and calibration. | Identified/ambiguous peaks, calibrated amounts/composition, recovery and detection limits. Fixed retention-time dictionaries and invented calibration factors cannot identify or quantify a new run. |

## Flow cytometry

Use [FlowKit](https://flowkit.readthedocs.io/en/latest/notebooks/flowkit-tutorial-part04-gates-module.html)
with the actual FCS channel metadata, compensation matrix, transformations and
gate hierarchy. Inspect acquisition-time, debris, singlet and viability gates
before phenotype gates; establish positivity from appropriate controls. Save
gates, per-sample event counts and both parent/total-population percentages.

- **Immunophenotyping / CD4 cytokines / senescence and apoptosis:** match the
  marker to the channel explicitly. Annexin/viability quadrants and SA-β-gal
  positivity need controls. Do not infer channel identity from the first
  fluorescence channel or define positives as the top 10/20% of every sample.
  Subtract unstimulated cytokine background only with matching design and
  denominators; retain raw frequencies too.
- **FACS selection:** selecting events by a gate is computational filtering,
  not physical cell sorting. Report selected event IDs, counts and fractions.
- **CFSE proliferation:** use a measured undivided control and a justified
  generation model. Define division/proliferation indices explicitly; recover
  precursor counts by the model's generation weighting before calculating
  precursor-based indices. A failed peak fit must not return the upstream
  fabricated division index 1.2 or proliferating fraction 65%.
- **Phase-duration estimation:** dual-pulse timing, labelling assumptions,
  controls and time-course counts are required for an identified kinetic model.
  Report fit uncertainty/sensitivity. One distribution of labels does not
  uniquely identify G1/S/G2-M durations or cell death.

Native fitting reference: [SciPy curve_fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html).
Choose an error model/weights appropriate to the measurements, and use sample-level
replicates rather than treating wells, frames or events as independent organisms.

Migration source: Biomni commit `400c1f366b96a35ca253e13c9b06c5076af41d65`,
assay functions in `biochemistry.py`, `cell_biology.py`, `immunology.py`,
`microbiology.py`, `pathology.py`, `physiology.py`, and `synthetic_biology.py`.
The old `isolate_purify_immune_cells` generates random yield/purity/viability;
it is retired. Protocol planning uses `lab-protocols`; measured outcomes use
the methods above.

Run fitting and statistics scripts with `biov exec python -- analysis.py`.
Run FlowKit scripts with `biov exec flowkit -- analysis.py`.

Read and write files in the working directory or configured data directories.
For a specialized tool without BioV, first obtain input files with BioV in a
`biov run prepare.py` script in the BioV application environment, then pass their paths to the tool.
