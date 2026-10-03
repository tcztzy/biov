---
name: analyze-physiological-signals
description: Analyze measured physiological time series, calcium traces, neural decoding data, and diffusion MRI. Use for ABR peaks, ciliary frequency, cosinor rhythms, hemodynamics, calcium events, or ADC maps.
---

# Physiological signal analysis

Require signal units, sampling times/rate, acquisition metadata, intervention
times and replicate IDs. Check missing samples, saturation and synchronization.
Keep raw and processed traces, processing parameters and exclusions. Use SciPy,
statsmodels, scikit-learn or the relevant imaging library directly; save numeric
results as TSV/Parquet/NIfTI plus a JSON summary and diagnostic plots.

| Task | Method | Outputs and restrictions |
| --- | --- | --- |
| ABR P1 | Use a predeclared species/acquisition-specific latency window and baseline convention; detect candidate peaks and inspect reproducibility across traces. | P1 latency, amplitude relative to stated baseline/trough and uncertainty. Do not select the largest peak anywhere or impose the upstream fixed window on a different acquisition. |
| Ciliary beat | Inspect the actual time-series ROI/video timestamps, detrend/window the signal and use a PSD/FFT method with sufficient duration. | Peak frequency, spectral width/power, temporal resolution and rejected ROIs. Respect Nyquist and frame timing; equally spaced arbitrary image tiles are not verified cilia ROIs. |
| Cosinor | Fit mesor plus sine/cosine terms at the prespecified period, with subject effects if repeated subjects are present. | Amplitude, phase in a declared time convention, intervals, rhythm test and residuals. Phase is unstable when amplitude is near zero; R² alone does not establish rhythmicity. |
| Calcium movies / endolysosomal traces | Use validated ROIs and motion/background correction; define baseline F0 and report ΔF/F0, peak/AUC/decay using actual times. | Traces, baseline interval, event rules, amplitudes, frequencies and decay-fit uncertainty. Event spacing is in seconds, not an undocumented number of frames; fit enough post-peak samples. Luminescence changes are not absolute calcium concentration. |
| Rhod-2 calcium | Require a real indicator calibration including Kd, Fmin/Fmax and the measurement conditions; otherwise report corrected relative fluorescence. | Calibrated concentration only inside the calibration range. A control image is not automatically Fmin and an image maximum is not Fmax. |
| Arterial pressure | Preserve the calibrated DC pressure level, inspect beats with a filter suitable for peak detection, then measure pressure values from the appropriate calibrated signal. | SBP, DBP, integrated MAP, pulse intervals and HR with artifact exclusions. A band-pass filtered trace has lost its baseline and cannot directly provide absolute SBP/DBP. |
| Neural behavior decoding | Align neural activity to behavior and fit supervised encoding/decoding on training trials; fit PCA/normalization inside the training split. Use blocked/trial or subject held-out evaluation and a baseline. | Predictions, error by held-out trial and selected hyperparameters. Randomly splitting autocorrelated time samples leaks information; an unsupervised Kalman state is not automatically the measured behavior. |
| ADC maps | Supply 4D diffusion data, b-values, orientation, mask and preprocessing history. Correct motion/distortion as required and fit the chosen diffusion signal model (for example with DIPY). | ADC and fit/QC maps in the source physical grid, with diffusivity units and invalid voxels marked. Match each volume to its b-value; a failed voxel fit is not zero diffusivity. |

Relevant APIs: [SciPy signal processing](https://docs.scipy.org/doc/scipy/reference/signal.html),
[statsmodels regression](https://www.statsmodels.org/stable/regression.html),
[scikit-learn model evaluation](https://scikit-learn.org/stable/modules/cross_validation.html),
and [DIPY reconstruction](https://docs.dipy.org/stable/examples_built/reconstruction/index.html).
These instructions support research analysis; they do not establish diagnostic
performance or reproduce the old function signatures.

Migration source: Biomni commit `400c1f366b96a35ca253e13c9b06c5076af41d65`,
signal functions in `physiology.py`, `bioengineering.py`, and `pathology.py`.

Run signal, statistical and DIPY analyses with
`biov exec python -- analysis.py`.

Read and write files in the working directory or configured data directories.
