---
name: model-biological-systems
description: Build and analyze explicitly specified biological mathematical models with SciPy, COBRApy, libSBML, and msprime. Use for kinetic networks, metabolic fluxes, population dynamics, demographic simulations, or parameter sweeps.
---

# Biological modelling

Begin with the actual equations/model file, state variables, units, parameter
sources, initial conditions and experimental question. Distinguish calibrated
models from illustrative simulations. Use native solvers; save SBML/JSON model
definitions, parameter tables and numeric trajectories (CSV/Parquet/NPZ without
object arrays), with a JSON summary naming versions, solver tolerances and seeds.
Never invent biological rate constants or claim that an assumed model predicts
an unmeasured experiment.

## Deterministic and stochastic models

| Model | Execution | Required checks |
| --- | --- | --- |
| General ODE / cell circuit / signalling | Use `scipy.integrate.solve_ivp` with a supplied RHS and parameters, explicit time/state units, solver and tolerances. For Hill/logic models state activation/inhibition rules and normalized versus concentration states. | Solver success, time coverage, physical bounds and conservation where applicable. The upstream four-state toy model is not a whole-cell model. Store equations, not an opaque callable/pickle. |
| Dimerization | Define homo/heterodimer stoichiometry and association/dissociation rates, or a properly constrained equilibrium model with Kd and total abundances. | Conservation of each monomer's total mass. An equilibrium affinity does not uniquely specify both kinetic rate constants. |
| Metabolic perturbation | Require kinetic laws and rate constants for a kinetic ODE model; represent interventions as explicit events/piecewise intervals. | State continuity at interventions, concentration units and mass balance. An SBML constraint model with arbitrary unit rate constants is not a validated kinetic model. Do not mutate state/parameters inside an adaptive solver's RHS. |
| RAS / thyroid compartments | Use supplied, referenced compartment equations, volumes, transport, binding and clearance parameters. Convert amounts/concentrations consistently across compartments. | Mass balance, parameter identifiability and sensitivity. Toy equations are educational simulations, not human dosing or physiological prediction. |
| Bacterial growth/clearance / gLV | Supply initial populations, growth and interaction/clearance parameters and carrying-capacity convention. Integrate the stated system; compare multiple initial conditions where stability matters. | Nonnegative states, integration success, equilibria and sensitivity. Fitting parameters from measured curves is a separate inference step. |
| Stochastic microbial populations | Use a maintained reaction-network simulator such as GillesPy2 with explicit birth/death propensities and integer initial counts. Run independently seeded replicates. | Propensities must stay nonnegative, including above carrying capacity. Report extinction counts and denominators separately from the times conditional on extinction. Label confidence intervals and Monte Carlo uncertainty correctly. |
| Demographic history | Construct `msprime.Demography`, validate event order, then call `sim_ancestry` and `sim_mutations` with explicit population sizes, ploidy/sample interpretation, rates, model and seeds. | Tree sequence plus VCF, sample IDs and parameter JSON. Event times are generations into the past; a Beta coalescent has its own model parameter assumptions. |

## SBML and flux balance

Use libSBML to create models with compartments, species, units, stoichiometry,
boundary conditions and supported kinetic laws. Parse formulas with libSBML,
run consistency checks and resolve model errors before exporting. SBML files
are the reusable interface; do not add a second BioV model representation.

For FBA, load a documented COBRApy SBML/JSON model; set exchange bounds,
medium, constraints and the actual objective reaction. Check solver status
before reading objective/fluxes. Report units and run flux variability or
knockout comparisons when needed by the question. A single optimal solution
does not establish a unique internal flux or a kinetic response.

For bifurcation questions, sweep the named parameter using declared initial
conditions and transients, solve equilibria/periodic states and assess their
stability with a suitable continuation method. A peak-count plot is descriptive;
mean log differences between adjacent samples is not a Lyapunov exponent and
must not establish chaos.

For digestion/process optimization, fit or import a mechanistic/empirical model
from real observations and bound the optimization to its supported domain.
Report validation and prediction uncertainty. Retire the upstream arbitrary
response surface; its optimum is an optimum of invented coefficients.

Official interfaces: [SciPy solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html),
[COBRApy](https://cobrapy.readthedocs.io/),
[libSBML](https://sbml.org/software/libsbml/),
[msprime](https://tskit.dev/msprime/docs/stable/), and
[GillesPy2](https://gillespy2.readthedocs.io/).

Migration source: Biomni commit `400c1f366b96a35ca253e13c9b06c5076af41d65`,
modelling functions in `systems_biology.py`, `synthetic_biology.py`,
`microbiology.py`, `genetics.py`, `bioengineering.py`, and `physiology.py`.
The skill preserves modelling tasks, not the upstream hard-coded equations
as validated biological models.

Run SciPy, libSBML, msprime and GillesPy2 scripts with
`biov exec python -- analysis.py`. Run COBRApy scripts with
`biov exec cobra -- flux_analysis.py`.

Read and write files in the working directory or configured data directories.
For a specialized tool without BioV, first obtain input files with BioV in a
`biov run prepare.py` script in the BioV application environment, then pass their paths to the tool.
