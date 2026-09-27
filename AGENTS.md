# decontX: notes for coding agents

The shared development standards (r-bioc-dev-standards) load automatically
at session start. This file adds only what is specific to this package.

## About

decontX implements two Bayesian decontamination methods for single-cell
data: **DecontX** (Yang et al. 2020) estimates and removes ambient RNA
contamination in scRNA-seq without needing empty droplets, and
**DecontPro** (Yin et al. 2023) removes ambient and background
contamination from CITE-seq ADT (protein) counts. Distributed through
Bioconductor since 3.18. There is no Shiny app and no pkgdown site yet
(see `dev/ROADMAP.md`), so those parts of the standards don't apply,
including `make site-check`.

## Layout

- `R/decon.R`: DecontX core. The `decontX()` S4 generic and methods, the
  EM algorithm, the `decontXcounts()` accessors, `simulateContamination()`,
  `retrieveFeatureIndex()`.
- `R/decontPro.R`: the `decontPro()` S4 generic and methods (Stan-based).
- `R/plot_decontx.R`, `R/plot_decontPro.R`: plotting functions.
- `R/celda_functions.R`: internal helpers copied from celda (see Related
  packages).
- `R/stan_helpers.R`: Stan model wiring (`.call_stan_vb()`,
  `.process_stan_vb_out()`).
- `src/DecontX.cpp` (Rcpp/RcppEigen EM steps) and `src/matrixSums*`
  (fast matrix operations).
- `inst/stan/shrinkage.stan`: the DecontPro Stan model. Changing it is a
  structural change and needs an ADR first.
- Generated files, never hand-edited:
  - Rcpp: `R/RcppExports.R`, `src/RcppExports.cpp`.
  - rstantools: `R/stanmodels.R`, `src/stanExports_shrinkage.*` (from
    `inst/stan/shrinkage.stan`).
- `man/examples/decontX.R` is **hand-written**. It is the shared example
  pulled into several help pages with `@example`, so edit it directly.
- `vignettes/decontX.Rmd`, `vignettes/decontPro.Rmd`.
- Style: the code uses **2-space indentation**, set in `.lintr`
  (`indentation_linter(indent = 2L)`). With 4 spaces, lintr would report
  over 800 indentation lints, against 3 with 2. Converting to 4 spaces
  would be a separate maintainer decision and its own PR.

## Object model

No custom S4 classes. The generics `decontX()`, `decontPro()` and
`decontXcounts()`/`decontXcounts<-` dispatch on `SingleCellExperiment`,
`Seurat` (decontPro only) and `ANY` (plain or sparse matrices). DecontX
stores its results in the input SingleCellExperiment:
- decontaminated counts in the `decontXcounts` assay (read it with
  `decontXcounts()`, not `assays()$`)
- contamination in `colData(sce)$decontX_contamination`
- cluster labels in `colData(sce)$decontX_clusters`
- the UMAP in `reducedDim(sce, "decontX_UMAP")`

`decontX(legacyInit = TRUE)` reproduces the pre-scrapper initialization.
It is a permanent backwards-compatibility option, kept for as long as the
upstream functions exist, so never describe it as transitional.

## Tests

- `tests/testthat/`: `test-decon.R`, `test-decontPro.R`,
  `test-matrixSums.R`, and shared fixtures in `helper-decontx.R`
  (testthat edition 3).
- Several decontPro tests run a real Stan variational fit. They use a
  10 x 8 fixture to stay under a second each; keep decontPro fixtures
  that small.
- Optional-package skips: the Seurat method test needs SeuratObject, and
  the `legacyInit` test needs scater.

## Extra make targets

The standard targets come from the shared `standards.mk` in
r-bioc-dev-standards; the Makefile holds only these extras. None are
allow-listed, so each one prompts.

- `coverage-report`: HTML coverage report (`covr::report()`); opens a
  browser. Ask first.
- `build`: builds the source tarball in the repo root. Ask first; the
  standard targets build in a temporary directory instead.
- `clean`: deletes compiled objects and tarballs. People only; denied in
  `.claude/settings.json`.

## Setup in a new worktree

- Install every package in Suggests (`devtools::install_dev_deps()`);
  checks don't relax missing Suggests.
- The first build compiles the Stan model and the Rcpp code, which takes
  several minutes. Later `devtools::load_all()` calls reuse the objects.
- `devtools::load_all()` and `make test` regenerate
  `src/stanExports_shrinkage.*`. Restore them (`git restore src/`) before
  committing.
- `make check` needs the `checkbashisms` script to check the rstantools
  `configure` scripts (`brew install checkbashisms` on macOS). Without it,
  R CMD check emits a WARNING.

## Related packages

- **celda:** `R/celda_functions.R` copies helpers from celda; keep the two
  copies from drifting and coordinate cross-package changes through
  issues. celda's `R/decon.R` is now a thin layer that forwards to
  `decontX::decontX()` (decontX is in celda's Suggests), so changes to
  decontX's exported API or results affect celda.
- **singleCellTK** imports celda, so it reaches DecontX through celda.
- ADR-0003 (`dev/adr/`) plans to reverse the celda/decontX dependency.

## Overrides

None.
