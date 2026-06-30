# CellSweeper 0.99.0

* Initial Bioconductor submission.
* Three-level QC framework for single-cell resolution spatial transcriptomics:
    * Level 1: Global pre-filtering (`globalFilter()`).
    * Level 2: Cluster-level artifact detection (`clusterCellTypes()`,
      `flagArtifactClusters()`).
    * Level 3: Within-cluster spatial QC (`clusterLocalOutliers()`).
* Convenience wrapper (`runCellSweeper()`).
* Visualization functions (`plotSpatialOutliers()`, `plotClusterQC()`,
  `plotSpatialDispersion()`, `plotClusterSummary()`).

## Pre-submission polish (0.99.0)

* `clusterLocalOutliers()` modified-z-score now offers two well-defined
  conventions via a new `robust_z` argument: `"mad"` (default,
  `(x - median(x)) / mad(x)` with scaled MAD) and `"iglewicz"`
  (`0.6745 * (x - median(x)) / mad(x, constant = 1)`). The previous
  default mixed the prefactor with the scaled MAD and was ~1.48x too
  conservative relative to the documented threshold.
* `permuteHomogeneity()` is now reproducible under any `BPPARAM`, using
  `BiocParallel::bpRNGseed()` to seed workers.
* `flagArtifactClusters(require_both = TRUE)` no longer silently
  disables artifact detection for samples with fewer than three
  clusters (where Mahalanobis cannot be computed); it falls back to
  the spatial dispersion flag alone for those samples.
* `clusterCellTypes()`:
  - Covariate regression is now a single vectorised `lm.fit()` call
    instead of a gene-by-gene loop, sharply reducing memory pressure
    on Xenium / Visium HD-scale data.
  - Now accepts any `SingleCellExperiment` (was: SPE-only).
  - Default `BSPARAM` is now `BiocSingular::IrlbaParam()` for tractable
    PCA cost at single-cell-resolution scale; PCA call uses the
    current `assay.type` argument.
  - RNG seeds for Harmony and the sketch-sampling step are scoped via
    `withr::with_seed()` so the user's RNG state is preserved.
* `runCellSweeper()` now forwards extra arguments to sub-functions via
  `...` instead of duplicating every default; `match.arg()` is used
  consistently across stages.
* `inst/CITATION` adds Stephanie Hicks and Boyi Guo (also added to
  `Authors@R`) and uses a hardcoded version string.
* Vignette `getting_started.Rmd` is now executed end-to-end against a
  small synthetic SpatialExperiment.
