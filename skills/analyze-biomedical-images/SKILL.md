---
name: analyze-biomedical-images
description: Segment, register, track, and quantify microscopy or biomedical images using established imaging tools. Use for cell/tissue morphology, multiplex imaging, histology, motility, colocalization, MRI registration, or volumetric measurements.
---

# Biomedical image analysis

Inspect axes, channel labels, bit depth, voxel spacing, timestamps and image
orientation before analysis. Require the specimen/ROI and intended measurement.
Keep original intensities for quantification: a contrast-enhanced display image
is not an intensity measurement. Save labelled masks, transforms or tracks plus
TSV/Parquet measurements and JSON parameters/QC; retain physical units.

## Segmentation and measurements

Use scikit-image for deterministic segmentation/measurement and Cellpose or
nnU-Net when a compatible trained model is justified. Check representative
masks against the source image before aggregating results. Record the model
and weight version, segmentation parameters, exclusions and biological sample ID.
An empty mask or failed load is a failure/empty result, never a synthetic image.

| Task | Method and required input | Measurements and limits |
| --- | --- | --- |
| Cells / microbial cells / colonies | Explicit channel, object size, ROI and pixel scale; threshold/seeded watershed or an appropriate Cellpose model. | Per-object area, axes, perimeter and intensity using `regionprops_table`; count touching objects only if the segmentation separates them. Colony count alone does not give CFU/mL without plated volume and dilution. |
| Myofibers / multiplexed tissue | Named nuclear and membrane/cytoplasmic channels; validate channel order and cell boundaries. | Cell/fiber masks, area and per-cell marker measurements. Expanded nuclear masks approximate cell boundaries; retain that limitation and excluded edge cells. |
| Cytoskeleton / mitochondrial morphology | Segment the specific structural marker, then measure object shape and skeleton topology using scikit-image/Skan. | Branch length, nodes and shape with calibration; connected components are not skeleton branches. Membrane-potential dye intensity is relative signal unless calibrated with suitable controls. |
| Corneal fibers | Marker-specific mask, physical scale and a stated 2D/3D ROI. | Area fraction or calibrated skeleton length in 2D; volume only from volumetric data. Do not label area fraction as nerve volume. |
| Plaques / CNS lesions / IHC / thrombus histology | Preserve stain channels, use stain separation where appropriate, and validate labels against annotated regions. | Counts, areas, intensity and spatial measurements. Arbitrary RGB cutoffs cannot establish lesion severity, cell identity, thrombus age or a clinical score. |
| Cell-cycle morphology | Extract measured shape/septum features; classification requires an organism-specific validated classifier or marker assay. | Feature table and, if justified, labels with validation. Size and intensity percentiles alone do not identify G1/S/G2/M. |
| Aortic geometry | Use a validated vessel/lumen mask and acquisition geometry. Measure cross-sections normal to the centreline when appropriate. | Calibrated diameters/areas and landmark definitions. The largest bright contour is not necessarily the aorta. |
| Bone micro-CT | Require true 3D acquisition, voxel spacing, tissue ROI and justified bone threshold; use a validated morphometry method such as BoneJ for local thickness. | BV/TV, thickness/separation and method-specific topology measures. Do not turn a 2D slice into 3D by adding an axis or call twice the mean distance-transform value validated trabecular thickness. |
| Blots / gels | Define lanes, band/background ROIs, detector linear range and loading controls on raw images. | Background-subtracted integrated intensities and normalization by sample. Histograms/automatic ROI suggestions need inspection; saturated bands cannot support quantitative ratios. |
| Colocalization | Register channels, exclude background and define the ROI/masks; use documented Pearson/Manders implementations. | Both coefficients with threshold/mask rules and biological replicate summaries. Intensity rescaling or selecting only double-positive pixels changes the estimand; colocalization does not prove molecular binding. |

## Motion and deformation

- **Cell migration and immune-cell flow:** detect objects, link with Trackpy using
  a search range based on displacement per frame, inspect tracks for swaps and
  gaps, then compute path length, net displacement, speed and persistence in
  physical units. Use actual frame/time differences, including missed frames.
  Mean squared displacement is the mean of squared displacement at each lag,
  not mean distance from the first point. Cluster standardized track features
  only when sample size and feature variance support it. Flow-relative rolling
  or arrest needs experiment-specific speed/duration definitions, not the
  upstream hard-coded five-frame threshold.
- **Tissue optical flow:** use OpenCV/scikit-image with image registration,
  physical spacing and time intervals. Retain the displacement vector field and
  valid-point mask. Calculate spatial derivatives on a spatial grid; reshaping
  sparse feature displacements into an arbitrary square array is invalid.
  State whether reporting velocity gradients, infinitesimal strain or finite
  deformation, and validate against known displacement.
- **Calcium movies:** motion-correct and define cell ROIs before extracting
  fluorescence. Use `analyze-physiological-signals` for baseline, event and decay
  calculations. Segmentation does not supply the missing sampling interval.

## Medical volumes and registration

1. Read NIfTI with NiBabel or SimpleITK and preserve affine/origin, spacing and
   direction. Split modalities only after identifying the fourth axis; it may
   represent time rather than the four BraTS contrasts assumed upstream.
2. For nnU-Net, use its native data naming/channel conventions, the selected
   dataset/configuration/folds, and model-specific preprocessing. Use verified
   weights and their usage terms; do not automatically download an unrelated
   model or globally patch `torch.load` to enable unsafe deserialization.
3. For registration use SimpleITK's rigid, affine or deformable registration
   with an appropriate similarity metric, optimizer and multiresolution setup.
   Preserve fixed/moving direction and transform files. Use nearest-neighbour
   interpolation for labels and suitable continuous interpolation for images.
4. Inspect aligned landmarks, checkerboards/overlays and deformable Jacobians
   where relevant; intensity correlation alone is not registration accuracy.
   Batch processing keeps a per-image transform and QC result and fails on a
   failed image instead of quietly switching transform types.
5. For segmentation overlays or surface extraction, use labels on the correct
   physical grid and marching cubes with voxel spacing plus affine transformation.
   A threshold-derived MRI surface is not validated facial anatomy.

Official methods: [scikit-image measurements](https://scikit-image.org/docs/stable/api/skimage.measure.html),
[Trackpy](https://soft-matter.github.io/trackpy/),
[SimpleITK registration](https://simpleitk.readthedocs.io/en/master/registrationOverview.html),
[nnU-Net](https://github.com/MIC-DKFZ/nnUNet),
[Cellpose](https://cellpose.readthedocs.io/), and [BoneJ](https://bonej.org/).
Validate the selected model on the user's imaging domain before interpreting
its measurements.

Migration source: Biomni commit `400c1f366b96a35ca253e13c9b06c5076af41d65`,
`bioimaging.py` and image-analysis functions in `bioengineering.py`,
`biophysics.py`, `cell_biology.py`, `immunology.py`, `microbiology.py`,
`pathology.py`, `pharmacology.py`, and `physiology.py`.

Run measurement, tracking and registration scripts with
`biov exec python -- analysis.py`. For the selected segmentation model, use
`biov exec nnunet -- ...` or `biov exec cellpose -- ...` with its native
prediction arguments.
nnU-Net uses v1 models and prediction options. Cellpose uses its current model
format; old cyto or Omnipose model names do not establish model compatibility.

Read and write files in the working directory or configured data directories.
For a specialized tool without BioV, first obtain input files with BioV in a
`biov run prepare.py` script in the BioV application environment, then pass their paths to the tool.
