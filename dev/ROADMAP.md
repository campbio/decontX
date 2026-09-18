# decontX Roadmap

Where the package is going. Proposals (from humans or agents) should
lean in these directions; anything structural still needs an ADR first.

## Near term

- **Quality baseline**: get `make bioccheck` clean and flip the CI
  BiocCheck job to a required check; raise test coverage, especially
  around `decontPro()` and the plotting functions.
- **Triage backlog**: convert existing BiocCheck/lintr findings into
  GitHub issues and burn them down incrementally (no big-bang cleanup).

## Medium term

- **pkgdown site**: decontX currently has no pkgdown site. Add
  `_pkgdown.yml`, build locally, deploy to `gh-pages` per the lab
  playbook (CI verifies structure with `pkgdown::check_pkgdown()` only —
  it never builds the site).
- **Precomputed vignettes**: the vignettes download data (TENxPBMCData,
  SingleCellMultiModal) and run MCMC — adopt the `.Rmd.orig` →
  precomputed `.Rmd` pattern so builds are fast and deterministic
  (record as an ADR when done).

## Longer term / ideas

- Performance work on the DecontX EM loop (C++/RcppParallel) guided by
  profiling, not speculation.
- Better interoperability helpers (e.g., round-tripping results through
  Seurat and anndata-based workflows).
- Keep `R/celda_functions.R` in sync with the celda package, or factor
  shared helpers into a common dependency (ADR territory).

## Non-goals (for now)

- No Shiny app in this package (interactive use is served by
  singleCellTK).
- No change to the statistical models without a methods-level review —
  model changes are science, not refactoring.
