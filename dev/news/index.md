# Changelog

## decontX 1.11.38

- A new “Getting started” vignette,
  [`vignette("decontX", package = "decontX")`](https://campbio.github.io/decontX/dev/articles/decontX.md),
  gives a short overview of DecontX and DecontPro, when to use each, and
  a quick example on simulated data. The DecontX workflow vignette that
  had that name is now
  [`vignette("decontX_RNA", package = "decontX")`](https://campbio.github.io/decontX/dev/articles/decontX_RNA.md);
  its workflow is unchanged
  ([\#84](https://github.com/campbio/decontX/issues/84)).
- decontX has a documentation website,
  <https://campbio.github.io/decontX/>, for the current Bioconductor
  release, and <https://campbio.github.io/decontX/dev/> for the
  development version. Its address is in DESCRIPTION and the README
  ([\#84](https://github.com/campbio/decontX/issues/84)).

## decontX 1.11.37

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  on a SingleCellExperiment with several batches now names the columns
  of each `decontX_<batch>_UMAP` `DecontX_UMAP_1` and `DecontX_UMAP_2`,
  as in a single-batch run’s `decontX_UMAP`. They were `UMAP_1` and
  `UMAP_2`, so code that reads these columns by name needs updating, and
  [`plotDecontXContamination()`](https://campbio.github.io/decontX/dev/reference/plotDecontXContamination.md)
  on one of them now titles its axes `DecontX_UMAP_1` and
  `DecontX_UMAP_2`
  ([\#99](https://github.com/campbio/decontX/issues/99)).
- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  now names its estimated parameters by ADT and droplet, from the row
  and column names of `filtered_counts`: the rows and columns of
  `delta_est` and `background_est`, the columns of `r_est`, and
  `delta_mean_est` and `background_mean_est`, which are now plain named
  vectors instead of 1-d arrays. Their values are unchanged
  ([\#129](https://github.com/campbio/decontX/issues/129)).
- [`plotDensity()`](https://campbio.github.io/decontX/dev/reference/plotDensity.md)
  and
  [`plotBoxByCluster()`](https://campbio.github.io/decontX/dev/reference/plotBoxByCluster.md)
  now stop when `counts` and `decontaminated_counts` have different
  numbers of droplets, or column names that differ. They paired the two
  matrices by position, so they gave a wrong plot with no error
  ([\#125](https://github.com/campbio/decontX/issues/125)).
- The README installs the development version from the `devel` branch.
  `install_github("campbio/decontX")` installs the current release
  ([\#137](https://github.com/campbio/decontX/issues/137)).

## decontX 1.11.36

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  now matches the barcodes of `background` to the cells within each
  batch when `batch` and `bgBatch` are given. Before, a barcode that was
  a cell in one batch also removed the empty droplet with the same
  barcode in another, as happens with bare 10x barcodes from several
  samples. Results change slightly for such data
  ([\#134](https://github.com/campbio/decontX/issues/134)).
- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  and
  [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  now warn when the called cells can’t be found in `background` by
  barcode: when the counts have no column names (only `background`
  without them warned before), and when no barcode matches but
  `background` has many droplets with as many counts as a typical cell
  (at least a tenth as many as there are cells, at or above the median
  cell’s total), as a raw matrix does when its barcodes are named
  differently from the cells’, for example without a sample prefix.
  These cells used to stay in the ambient profile without a word
  ([\#134](https://github.com/campbio/decontX/issues/134)).
- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  help and vignette say that `background` barcodes are matched across
  all droplets, so with several samples they must be unique across
  samples ([\#134](https://github.com/campbio/decontX/issues/134)).

## decontX 1.11.35

- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  now removes the called cells from `background`, as
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  does: droplets whose barcodes (column names) are also in
  `filtered_counts` are dropped, so the raw matrix can be given as it
  is. Before, their native signal went into the ambient profile, so
  results change for calls whose `background` included called cells. A
  `background` without column names gives a warning, and one with no
  droplets left, or no counts in the ADTs of `filtered_counts`, stops
  with a clear error instead of Stan’s “not a simplex” error
  ([\#131](https://github.com/campbio/decontX/issues/131)).

## decontX 1.11.34

- The decontPro vignette is rewritten. It uses the cord blood data’s own
  cell-type labels (no Seurat clustering), explains every output, shows
  how to tune the priors with before/after density plots, has a template
  for estimating the ambient profile from empty droplets with
  `background`, and covers run time, reproducibility, a FAQ and how to
  cite DecontPro. Stan’s progress no longer fills the page
  ([\#83](https://github.com/campbio/decontX/issues/83)).
- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  help describes every element it returns and says that only the ADT
  totals of `background` are used, so it must hold empty droplets only.
  Its example is smaller and seeded
  ([\#83](https://github.com/campbio/decontX/issues/83)).
- dplyr is no longer suggested
  ([\#83](https://github.com/campbio/decontX/issues/83)).

## decontX 1.11.33

- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  arguments are now camelCase, like decontX’s: `cellType`, `deltaSd`,
  `backgroundSd`, and `background` for the empty-droplet counts. The old
  names `cell_type`, `delta_sd`, `background_sd` and `ambient_counts`
  keep working as permanent aliases; giving both names of an argument is
  an error. Positional calls are unchanged. Two kinds of call that used
  to work now stop with an error: an abbreviated old name (such as
  `ambient =`, or `background =` for `background_sd`, since `background`
  now names the empty-droplet counts), and a positional value after an
  alias (such as `cell_type = k, 1e-2`)
  ([\#81](https://github.com/campbio/decontX/issues/81)).
- New
  [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  arguments `iter` (default 50000) and `seed` (default

12345. control Stan’s fit; the defaults are the values used before, so
       results don’t change. `verbose = FALSE` hides decontPro’s
       messages and Stan’s console output, which is shown only if the
       fit fails ([\#81](https://github.com/campbio/decontX/issues/81)).

- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  deletes the temporary CSV that rstan writes, about 100 bytes per ADT x
  droplet, as soon as the fit ends; before, it stayed until R exited
  ([\#46](https://github.com/campbio/decontX/issues/46),
  [\#81](https://github.com/campbio/decontX/issues/81)).
- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  help: `deltaSd` and `backgroundSd` are standard deviations of the
  priors, not variances, and `background` is described apart from the
  model’s background component
  ([\#81](https://github.com/campbio/decontX/issues/81)).

## decontX 1.11.32

- The decontX vignette shows the code that maps clusters to cell types,
  reuses the clusters for its second run, and splits “Running decontX”
  into runs without and with the background matrix, showing how to
  import the raw matrix and what decontX needs from it
  ([\#43](https://github.com/campbio/decontX/issues/43)). New sections
  cover multiple samples
  ([\#45](https://github.com/campbio/decontX/issues/45)) and a FAQ on
  data with one cell type
  ([\#43](https://github.com/campbio/decontX/issues/43)) and doublet
  detection ([\#44](https://github.com/campbio/decontX/issues/44)). It
  also says how to cite DecontX.
- The package help page, the README and DESCRIPTION cite both papers.
- Documentation fixes: contamination is a proportion between 0 and 1,
  not a percentage;
  [`simulateContamination()`](https://campbio.github.io/decontX/dev/reference/simulateContamination.md)
  gives its real `delta` default (`c(1, 10)`) and lists everything it
  returns;
  [`decontXcounts()`](https://campbio.github.io/decontX/dev/reference/decontXcounts.md)
  has an example ([\#83](https://github.com/campbio/decontX/issues/83)).

## decontX 1.11.31

- [`plotDensity()`](https://campbio.github.io/decontX/dev/reference/plotDensity.md)
  and
  [`plotBoxByCluster()`](https://campbio.github.io/decontX/dev/reference/plotBoxByCluster.md)
  now use a factor in `features` by its labels; before, they plotted the
  row at the factor’s integer code under the label’s title
  ([\#122](https://github.com/campbio/decontX/issues/122)). `features`
  must be rownames: row numbers, which used to select rows by position,
  now stop with an error. Empty `features`, features that aren’t
  rownames of both count matrices, and (in
  [`plotBoxByCluster()`](https://campbio.github.io/decontX/dev/reference/plotBoxByCluster.md))
  a `cell_type` without one entry per droplet now stop with an error
  that names the argument. Before, a `cell_type` of the wrong length
  could be recycled silently, giving droplets the wrong cell types
  ([\#122](https://github.com/campbio/decontX/issues/122),
  [\#123](https://github.com/campbio/decontX/issues/123)).

## decontX 1.11.29

- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  splits the counts after the Stan fit with about half the memory and in
  less time. The split results can differ from before in the last digits
  (about 1e-16, relative); the estimated parameters are unchanged
  ([\#111](https://github.com/campbio/decontX/issues/111)).

## decontX 1.11.28

- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  now stops with a clear error when the counts have a single ADT or a
  single droplet. Before, Stan failed with an unclear message about
  mismatched dimensions
  ([\#110](https://github.com/campbio/decontX/issues/110)).

## decontX 1.11.26

- decontX now needs scrapper 1.5.0 or later; every Bioconductor release
  since 3.23 ships one. Results are unchanged
  ([\#115](https://github.com/campbio/decontX/issues/115)).

## decontX 1.11.25

- Removed compiled code that nothing used, and the help pages of
  internal functions. The results of
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md),
  [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  and
  [`simulateContamination()`](https://campbio.github.io/decontX/dev/reference/simulateContamination.md)
  are unchanged ([\#78](https://github.com/campbio/decontX/issues/78)).

## decontX 1.11.24

- decontX no longer imports plyr or reshape2. Their few uses are now
  done in base R, with the same results
  ([\#77](https://github.com/campbio/decontX/issues/77)).

## decontX 1.11.23

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  now stops with an error when `batch` or `bgBatch` contains an empty
  string (`""`). The estimates of a batch named `""` couldn’t be read
  back by name, and its cells lost their cluster labels in `z` when
  decontX estimated the clusters or when the labels also appeared in
  another batch: every cell got `"-"` (or `""` when it was the only
  batch). With more than one batch, a SingleCellExperiment and the
  default `embedding = TRUE`, the run failed with an unclear error while
  storing the UMAPs. Contamination and decontaminated counts were right
  ([\#102](https://github.com/campbio/decontX/issues/102)).

## decontX 1.11.21

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  builds the decontaminated counts faster and with less memory: on
  100,000 simulated cells with `z` given and `embedding = FALSE`, the
  run takes 1.2 s instead of 2.8 s and peaks about 1.2 GB lower. The
  native probability is now computed in the same direct form as in the
  EM, so the decontaminated counts can change in the last digits (at
  most about 1e-14 relative); contamination estimates are unchanged
  ([\#80](https://github.com/campbio/decontX/issues/80)).

## decontX 1.11.20

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)’s
  EM step is faster: a run with cluster labels given in `z` and
  `embedding = FALSE` takes about half the time on the PBMC4k data (3.4
  s to 1.7 s). The results are unchanged
  ([\#80](https://github.com/campbio/decontX/issues/80)).

## decontX 1.11.19

- New
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  argument `embedding`. With `embedding = FALSE` and cluster labels
  given in `z`,
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  skips the embedding (the normalization, variable genes, PCA and UMAP),
  which is most of the run time on large data sets. The decontamination
  is the same (except with `legacyInit = TRUE` and `seed = NULL`, where
  the EM’s random start differs); only the decontX UMAP is left out of
  the results. `embedding = FALSE` needs `z`, and lowers the minimum
  batch size from 3 cells (15 with `legacyInit = TRUE`) to 2
  ([\#79](https://github.com/campbio/decontX/issues/79)).
- New
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  argument `numThreads` (default 1) for computing the embedding. With
  more than one thread the PCA can differ in the last digits on large
  data sets, which can move the UMAP coordinates and change estimated
  clusters: on 100,000 simulated cells, four threads found 11 clusters
  instead of 10. Set `numThreads` the same way to reproduce a run. With
  `legacyInit = TRUE` only the UMAP uses the threads
  ([\#79](https://github.com/campbio/decontX/issues/79)).
- [`plotDecontXContamination()`](https://campbio.github.io/decontX/dev/reference/plotDecontXContamination.md)
  can plot the contamination on any reducedDim of a SingleCellExperiment
  with the new argument `reducedDimName`, for example after a run with
  `embedding = FALSE`. With `reducedDimName` and `batch = NULL` it plots
  every cell ([\#79](https://github.com/campbio/decontX/issues/79)).

## decontX 1.11.17

- With HDF5-backed counts, `decontXcounts` is now a DelayedMatrix held
  in memory instead of one written to a temporary HDF5 file. Runs are
  faster, and a saved object no longer points at a file that is deleted
  when the R session ends. The result needs memory for the sparse
  decontaminated counts, and HDF5Array’s dump settings (such as
  [`setHDF5DumpDir()`](https://rdrr.io/pkg/HDF5Array/man/dump-management.html))
  no longer apply to it. To keep it on disk, save it explicitly:
  [`HDF5Array::saveHDF5SummarizedExperiment()`](https://rdrr.io/pkg/HDF5Array/man/saveHDF5SummarizedExperiment.html)
  for a SingleCellExperiment, or
  [`HDF5Array::writeHDF5Array()`](https://rdrr.io/pkg/HDF5Array/man/writeHDF5Array.html)
  for the matrix. Values don’t change
  ([\#76](https://github.com/campbio/decontX/issues/76)).

## decontX 1.11.16

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  makes fewer copies of the decontaminated counts when it puts them
  together, so it allocates less memory and runs slightly faster, with
  or without `batch`. Results don’t change
  ([\#76](https://github.com/campbio/decontX/issues/76)).

## decontX 1.11.15

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  stops when `bgBatch` is given with `batch` but without `background`.
  It used to ignore `bgBatch` then, unless its values didn’t match
  `batch` ([\#76](https://github.com/campbio/decontX/issues/76)).

## decontX 1.11.14

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  on a matrix also returns `batch`, the batch of each cell used for the
  run (`"all_cells"` when `batch` isn’t given). The SingleCellExperiment
  results don’t change, apart from `runParams` below
  ([\#76](https://github.com/campbio/decontX/issues/76)).
- `runParams` in the
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  results always has `batchBackground`, which is `NULL` when
  `background` and `bgBatch` are both `NULL`. It used to be missing then
  ([\#96](https://github.com/campbio/decontX/issues/96)).

## decontX 1.11.13

- The marker plots’ error for duplicated `groupClusters` values has its
  missing space back (“assigned to a single group”)
  ([\#75](https://github.com/campbio/decontX/issues/75)).

## decontX 1.11.12

- [`simulateContamination()`](https://campbio.github.io/decontX/dev/reference/simulateContamination.md)
  gives every cell `NRange[1]` counts when both ends of `NRange` are
  equal; it used to draw totals from 1 to that value. `K = 1` now stops
  with a clear message: the simulated contamination comes from the other
  cell populations, so it needs at least two. So does a `C` too small to
  give two populations (for example `C = 1`). Both used to fail with
  “Division by 0”
  ([\#74](https://github.com/campbio/decontX/issues/74)).
- [`retrieveFeatureIndex()`](https://campbio.github.io/decontX/dev/reference/retrieveFeatureIndex.md)
  never matches a feature that is `NA`, even when a row is named `NA`,
  and stops with a clear message when `by` isn’t a single string
  ([\#74](https://github.com/campbio/decontX/issues/74)).
- `runParams` in the
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  results also records `seed` and `iterLogLik`
  ([\#74](https://github.com/campbio/decontX/issues/74)).
- Typos fixed in messages and help pages. The error for `bgBatch` given
  without `batch` now names both arguments, and the
  [`simulateContamination()`](https://campbio.github.io/decontX/dev/reference/simulateContamination.md)
  help page names its outputs `nativeCounts` and `observedCounts`
  ([\#74](https://github.com/campbio/decontX/issues/74)).

## decontX 1.11.11

- [`plotDensity()`](https://campbio.github.io/decontX/dev/reference/plotDensity.md)
  smooths the original and decontaminated curves of an ADT with one
  bandwidth, computed from both. Each curve used to get its own, which
  made the less spread-out curve’s peaks look taller than they are
  ([\#73](https://github.com/campbio/decontX/issues/73)).
- [`plotDensity()`](https://campbio.github.io/decontX/dev/reference/plotDensity.md)
  and
  [`plotBoxByCluster()`](https://campbio.github.io/decontX/dev/reference/plotBoxByCluster.md)
  save to `file` as given when it ends in “.pdf”; `file = "x.pdf"` used
  to write `x.pdf.pdf`. Without the extension, “.pdf” is still added.
  The help pages now say the path is returned invisibly
  ([\#73](https://github.com/campbio/decontX/issues/73)).
- When `groupClusters` names labels that aren’t in `z`, the marker plots
  also list the labels `z` has. After a batched
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  run, those can carry the batch as a prefix, such as `"b1-1"`, which
  the help pages now explain
  ([\#73](https://github.com/campbio/decontX/issues/73)).
- [`plotDecontXContamination()`](https://campbio.github.io/decontX/dev/reference/plotDecontXContamination.md)
  reads the UMAP columns by position, so a UMAP without column names, as
  after an h5ad round trip, plots instead of failing with
  `.data[[NULL]]`
  ([\#73](https://github.com/campbio/decontX/issues/73)).
- The
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  help page says a re-run keeps UMAPs from an earlier run with different
  batches in `reducedDims`
  ([\#73](https://github.com/campbio/decontX/issues/73)).
- The decontPro vignette passes Seurat’s cluster identities to
  [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  as they are, so the cluster numbers in
  [`plotBoxByCluster()`](https://campbio.github.io/decontX/dev/reference/plotBoxByCluster.md)
  match the Seurat plots; they used to be shifted by one. Results don’t
  change ([\#73](https://github.com/campbio/decontX/issues/73)).

## decontX 1.11.10

- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  reads every counts layer of a Seurat v5 object whose assay is split,
  joining them in
  [`Cells()`](https://satijalab.github.io/seurat-object/reference/Cells.html)
  order. It used to read only the first layer, with a warning, and
  either stop in Stan or run on that layer’s cells alone
  ([\#72](https://github.com/campbio/decontX/issues/72)).
- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  accepts any `cell_type` labels, such as integers with gaps, character
  strings or a factor, and recodes them for Stan. The rows of
  `parameters$r_est` are named by label. A `cell_type` of the wrong
  length or with NA, counts that aren’t non-negative whole numbers, and
  NA or negative values in `ambient_counts` stop with a clear message,
  and a failed Stan fit says so instead of stopping with “no applicable
  method for `@`”. Results for labels 1 to K don’t change
  ([\#72](https://github.com/campbio/decontX/issues/72)).
- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  stops on arguments it doesn’t know, such as `seed`, `batch` or a
  misspelled `delta_SD`, which used to be ignored
  ([\#72](https://github.com/campbio/decontX/issues/72)).

## decontX 1.11.9

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  stops with a clear message when a batch, or `x` without `batch`, has
  fewer than 3 cells (15 with `legacyInit = TRUE`), which the
  initialization needs. A one-cell batch used to stop with “subscript
  out of bounds” ([\#71](https://github.com/campbio/decontX/issues/71)).
- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  accepts a DelayedMatrix over a dense matrix, which used to stop after
  the whole run, and dgRMatrix and dgTMatrix counts without errors or
  deprecation warnings. Other inputs, such as a data.frame, a
  SummarizedExperiment or a matrix of logical or character values
  (logical sparse matrices used to run as 0/1 counts), stop with a
  message naming the class. With a SingleCellExperiment `x`,
  `background` must be a SingleCellExperiment too
  ([\#71](https://github.com/campbio/decontX/issues/71)).
- `decontXcounts` is always a dgCMatrix, or a DelayedMatrix for
  DelayedArray input. dgTMatrix and SVT_SparseMatrix input used to come
  back in their own class. No results change
  ([\#71](https://github.com/campbio/decontX/issues/71)).

## decontX 1.11.8

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  returns `z` as a factor. A supplied `z` comes back with its own labels
  instead of integer codes; a factor keeps its level order. In batched
  runs the labels get the batch prefix (`"<batch>-<label>"`) only when a
  label is found in more than one batch. Estimated clusters keep their
  numbers (`"<batch>-<k>"` with `batch`). The `phi` and `eta` estimates
  are named by cluster label, and each batch’s `estimates` no longer
  holds `z` twice. `as.integer(result$z)` may give different values than
  before, and `groupClusters` in the marker plots now takes these labels
  for a supplied `z`
  ([\#70](https://github.com/campbio/decontX/issues/70)).
- Unused levels in a factor `z`, including levels missing from one
  batch, no longer give `NaN` columns in `phi`. Contamination for such a
  `z` changes by about 1e-16; no other results change
  ([\#70](https://github.com/campbio/decontX/issues/70)).
- [`plotDecontXMarkerPercentage()`](https://campbio.github.io/decontX/dev/reference/plotDecontXMarkerPercentage.md)
  gives the right percentages and labels for any cluster labels; labels
  other than 1 to K used to stop, give NA labels or mislabel clusters.
  Both marker plots work with one marker gene or a `groupClusters` that
  keeps one cell, and
  [`plotDecontXMarkerExpression()`](https://campbio.github.io/decontX/dev/reference/plotDecontXMarkerExpression.md)
  no longer assigns cells to the wrong cluster when the cell names look
  like numbers ([\#70](https://github.com/campbio/decontX/issues/70)).

## decontX 1.11.7

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  allows a `z` with one cluster when `background` is supplied, since the
  contamination profile then comes from the empty droplets. This fits
  only when the cells are a single cell type. Without `background`, one
  cluster still stops, and the error now names the batch
  ([\#69](https://github.com/campbio/decontX/issues/69)).
- When dbscan finds only one cluster, the k-means fallback now runs on
  the PCA components instead of a dense matrix of all genes, which could
  need gigabytes of memory
  ([\#43](https://github.com/campbio/decontX/issues/43)). This changes
  the clusters and contamination estimates on that path only;
  `legacyInit = TRUE` still clusters all genes.
- With `legacyInit = TRUE`, the k-means fallback no longer stops with
  “argument is not a matrix”. It had failed on this path since
  `legacyInit` was added, so no earlier results change.

## decontX 1.11.6

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  now matches the rows of `background` to the genes of `x` by name, and
  [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  matches the rows of `ambient_counts` to the ADTs of `filtered_counts`.
  Before, rows were matched by position: rows in a different order ran
  without an error and gave wrong estimates, extra or missing genes gave
  a cryptic error
  ([\#43](https://github.com/campbio/decontX/issues/43)), and names that
  didn’t overlap ran. Now extra rows are dropped with a log message,
  missing genes stop with an error that names them, and rows are matched
  by position only when either matrix has no row names. Results don’t
  change when the rows already match
  ([\#68](https://github.com/campbio/decontX/issues/68)).
- `bgBatch` must have one value per background droplet, with no `NA`,
  and its values must be the same as those of `batch`. Before, a short
  `bgBatch` was recycled, and a batch without matching background
  droplets silently ran without background.
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  also stops when no background droplets are left after removing
  barcodes that are in `x`, and a batch with a single background droplet
  no longer fails.

## decontX 1.11.5

- [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  now checks its inputs before running and stops with an error that
  names the argument. Before, a `z` or `batch` shorter than the number
  of cells was padded with `NA` and gave wrong results; `z` or `batch`
  with `NA`, invalid `maxIter`, `iterLogLik`, `convergence`, `seed`,
  `varGenes`, `dbscanEps` or `logfile` values, non-logical
  `estimateDelta`, `legacyInit` or `verbose`, a `delta` that isn’t two
  finite positive numbers, and negative counts in `x` or `background`
  ran silently or failed with cryptic errors. `maxIter`, `iterLogLik`
  and `varGenes` must be whole numbers (or `Inf`, which keeps working),
  and `maxIter = Inf` needs `convergence` above 0. An `assayName` or
  `bgAssayName` that isn’t an assay of the input now stops before
  decontX runs; an assay index still works. Results don’t change for
  valid input ([\#67](https://github.com/campbio/decontX/issues/67)).
- The `varGenes` and `legacyInit` docs now say they still set the UMAP
  stored in the results when `z` is given, and the error for cells with
  no counts no longer suggests supplying `z`.

## decontX 1.11.3

- The “Working with Seurat” section of the decontX vignette now runs as
  pasted with SeuratObject 5: it reads counts with
  [`LayerData()`](https://satijalab.github.io/seurat-object/reference/Layers.html)
  instead of `GetAssayData(slot = )` (defunct since SeuratObject 5.3),
  attaches SingleCellExperiment, and uses one name for the Seurat object
  ([\#65](https://github.com/campbio/decontX/issues/65)).

## decontX 1.11.2

- The `decontX-package` help page now has a Value section, which fixes a
  BiocCheck warning.
- Documentation is regenerated with roxygen2 8.1.0, and the Stan model
  exports with rstantools 2.7.1. Results don’t change.

## decontX 1.11.1

- `decontX` no longer depends on the `celda` package (ADR-0003 Stage 1).
  The one internal use of `celda::normalizeCounts()` (proportion
  normalization inside
  [`simulateContamination()`](https://campbio.github.io/decontX/dev/reference/simulateContamination.md))
  is replaced by the package’s own `fastNormProp` routine, which
  produces identical results, and the example and vignette now plot
  cluster labels and marker-gene expression on the UMAP with
  [`scater::plotReducedDim()`](https://rdrr.io/pkg/scater/man/plotReducedDim.html)
  instead of celda’s plotting functions. This drops `celda` (and its
  transitive dependencies) entirely and removes the `R CMD check`
  deprecation warnings that came from celda’s masked, stale copy of
  `decontX` being invoked after
  [`library(celda)`](https://rdrr.io/r/base/library.html) in the
  examples.
- The cell cluster initialization used by
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  when `z` is not supplied (log-normalization, variable gene selection,
  PCA, and UMAP) has been migrated from scater/scuttle to the scrapper
  package, because scuttle deprecated
  [`logNormCounts()`](https://rdrr.io/pkg/scuttle/man/logNormCounts.html)/[`normalizeCounts()`](https://rdrr.io/pkg/scrapper/man/normalizeCounts.html)
  in Bioconductor 3.24. This removes the deprecation warnings from
  `R CMD check`.
- Results of
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  without `z` may differ slightly from previous versions (same workflow
  and parameters; different variable-gene criterion and UMAP
  implementation). A new `legacyInit = TRUE` argument reproduces the
  previous behavior; it requires the scater package (now in Suggests)
  and is kept for backwards compatibility for as long as the upstream
  functions remain available. Users supplying their own `z` are
  unaffected.
- With the default initialization, the UMAP step itself is deterministic
  even when `seed = NULL` (scrapper seeds it explicitly), but overall
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  results with `seed = NULL` remain non-reproducible because the EM
  initialization is unseeded — set a seed for full reproducibility.
- The plotting functions
  ([`plotDecontXContamination()`](https://campbio.github.io/decontX/dev/reference/plotDecontXContamination.md),
  [`plotDecontXMarkerPercentage()`](https://campbio.github.io/decontX/dev/reference/plotDecontXMarkerPercentage.md),
  [`plotDecontXMarkerExpression()`](https://campbio.github.io/decontX/dev/reference/plotDecontXMarkerExpression.md),
  [`plotDensity()`](https://campbio.github.io/decontX/dev/reference/plotDensity.md),
  [`plotBoxByCluster()`](https://campbio.github.io/decontX/dev/reference/plotBoxByCluster.md))
  no longer use the deprecated
  [`ggplot2::aes_string()`](https://ggplot2.tidyverse.org/reference/aes_.html)
  or `trans =` scale argument, which removes the ggplot2 deprecation
  warnings they emitted. decontX now requires ggplot2 (\>= 3.5.0).
- [`plotDecontXMarkerPercentage()`](https://campbio.github.io/decontX/dev/reference/plotDecontXMarkerPercentage.md)
  and
  [`plotDecontXMarkerExpression()`](https://campbio.github.io/decontX/dev/reference/plotDecontXMarkerExpression.md)
  now work when `x` is a dense or sparse matrix (they previously failed
  with a ggplot2 aesthetics error when the plot was drawn). The facets
  of
  [`plotDecontXMarkerExpression()`](https://campbio.github.io/decontX/dev/reference/plotDecontXMarkerExpression.md)
  now follow the order of `groupClusters`, matching the x-axis.
- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  no longer draws 1000 unused samples from the variational posterior.
  Only the mean of the approximation is used, but the draws scaled with
  the number of ADTs x droplets and could take hundreds of GB of disk
  and memory on large datasets
  ([\#46](https://github.com/campbio/decontX/issues/46)). On a 50 x 400
  test dataset the temporary file shrank from 444 MB to 2 MB and peak
  memory from 3.0 GB to 1.2 GB. Results are unchanged. rstan’s “Pareto k
  diagnostic” warning is no longer shown: it is estimated from the draws
  and would always report `Inf` with the minimal number of draws now
  requested.

## decontX 1.10.0

- No significant changes (Bioconductor 3.23 release).

## decontX 1.8.0

- No significant changes (Bioconductor 3.22 release).

## decontX 1.6.0

- No significant changes (Bioconductor 3.21 release). Includes the
  DecontPro changes released in 1.4.1 on the Bioconductor 3.20 branch.

## decontX 1.4.1 (2025-02-16)

- [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  gained an `ambient_counts` argument, so an empty-droplet count matrix
  can be supplied to estimate the ambient contamination profile
  empirically instead of inferring it from the filtered droplets.
- **Breaking:** the first argument of
  [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  was renamed from `object` to `filtered_counts` for consistency with
  the new `ambient_counts` argument. Code that passed the count matrix
  by name (`decontPro(object = x, ...)`) must be updated; positional
  calls are unaffected.
- The estimated ambient contamination profile is now returned in the
  `parameters` element of the
  [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  output (`p_est`).
- Updated deprecated array syntax in the Stan model, fixed a C compiler
  warning in `matrixSums`, and corrected a `deconPro` typo in the
  DecontPro vignette.

## decontX 1.4.0

- No significant changes (Bioconductor 3.20 release).

## decontX 1.2.0

- No significant changes (Bioconductor 3.19 release).

## decontX 1.0.0

- First Bioconductor release (Bioconductor 3.18), providing
  [`decontX()`](https://campbio.github.io/decontX/dev/reference/decontX.md)
  for ambient RNA decontamination of single-cell RNA-seq counts and
  [`decontPro()`](https://campbio.github.io/decontX/dev/reference/decontPro.md)
  for single-cell protein expression (CITE-seq/TotalSeq) data.

## decontX 0.99.5 (2023-10-19)

- First submission to Bioconductor after review.
