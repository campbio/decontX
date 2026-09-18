# AGENTS.md

## Campbell Lab Playbook (common across lab packages — v2.0, do not edit per-repo)

### Common commands
`make test` / `make check` / `make bioccheck` / `make docs` / `make lint`
(See Makefile for definitions. These are the ONLY sanctioned entry points.
Run `make test` after every change; `make check` before opening a PR.)

### Git and PR workflow
- Branch from devel; all work lands via PR. Never push to devel or master.
- Use plan mode for any non-trivial change.
- Run /code-review before requesting human review.
- Every user-facing change gets a NEWS.md entry.

### Coding conventions
- Style enforced by lintr/styler (config in repo); <= 80-char lines (BiocCheck).
- roxygen2 owns man/ and NAMESPACE — NEVER hand-edit them.
- Use accessor functions, not @ slot access, outside class definition files.

### Documentation (pkgdown)
- The website is GENERATED. Improve docs by editing roxygen comments, vignettes,
  and _pkgdown.yml — never files under docs/ or the gh-pages branch.
- New exported functions MUST be added to the _pkgdown.yml reference index;
  verify with pkgdown::check_pkgdown().
- To preview one changed page: pkgdown::build_article("<name>") or
  build_reference_index(). NEVER run a full build_site() as verification —
  full site builds/deploys are a local maintainer action (make site-deploy).
- Files under vignettes/articles/ are pkgdown-only and NOT checked by
  R CMD check — knit locally when you edit them.

### Shiny app rules (packages with inst/shiny only)
- The app contains NO analysis logic. Server code only wires inputs to
  exported package functions and renders results. New app features are
  implemented as tested, exported functions first.
- Reactive logic is tested with shiny::testServer(); the golden path is
  covered by a small shinytest2 smoke suite (make test-app).
- UI changes are verified with a screenshot of the RUNNING app
  (make app + browser), not just passing tests.
- inst/ code is invisible to R CMD check — tests and lintr are the only
  guards; inst/shiny is included in the lint paths.

### Versioning and releases
- Bioconductor even/odd x.y.z scheme; releases ~April and ~October.
- Follow dev/RELEASE.md for the release checklist.

### Safety rules
- No structural refactors (file splits, DESCRIPTION dependency changes,
  class redesign) without an approved ADR — propose via a GitHub issue.
- Never commit secrets, tokens, or absolute local paths.
- Architectural decisions are recorded in dev/adr/ (see template and index
  there). Never store anything in docs/ — that is pkgdown build output.
- Maintainer docs (release, roadmap, audits) live in dev/, not the root.

## This package: decontX

### Project overview
decontX implements two Bayesian decontamination methods for single-cell
genomics data: **DecontX** (Yang et al. 2020), which estimates and removes
ambient RNA contamination in single-cell RNA-seq without requiring empty
droplet information, and **DecontPro** (Yin et al. 2023), which estimates
and removes contamination from ambient and background sources in CITE-seq
ADT (protein) data. It is distributed through Bioconductor (first release
in Bioc 3.18).

### Repository map
- `R/decon.R` — DecontX core: `decontX()` S4 generic + methods, the EM
  algorithm, `decontXcounts()` accessors, `simulateContamination()`,
  `retrieveFeatureIndex()`.
- `R/decontPro.R` — `decontPro()` S4 generic + methods (Stan-based).
- `R/plot_decontx.R`, `R/plot_decontPro.R` — plotting functions.
- `R/celda_functions.R` — internal helpers shared with the celda package.
- `R/stan_helpers.R`, `R/stanmodels.R` — Stan model wiring (rstantools).
- `src/` — C/C++ code: `DecontX.cpp` (Rcpp/RcppEigen EM steps),
  `matrixSums*` (fast matrix ops), `stanExports_shrinkage.*` (GENERATED
  from `inst/stan/shrinkage.stan` by rstantools).
- `inst/stan/shrinkage.stan` — the DecontPro Stan model source.
- `tests/testthat/` — `test-decon.R`, `test-decontPro.R`,
  `test-matrixSums.R` (testthat edition 3).
- `vignettes/decontX.Rmd`, `vignettes/decontPro.Rmd`.

### Object model
No custom S4 classes are defined. The package exposes S4 generics —
`decontX()`, `decontPro()`, `decontXcounts()`/`decontXcounts<-` — with
methods dispatching on `SingleCellExperiment`, `Seurat` (decontPro only),
and `ANY` (plain/sparse matrices). DecontX results are stored in the input
SingleCellExperiment: decontaminated counts in the `decontXcounts` assay
(use the `decontXcounts()` accessor, never `assays()$` directly),
contamination estimates in `colData(sce)$decontX_contamination`, cluster
labels in `colData(sce)$decontX_clusters`, and the UMAP in
`reducedDim(sce, "decontX_UMAP")`.

### Environment setup
- R >= 4.3.0 plus a working C++ toolchain and GNU make (SystemRequirements);
  macOS needs Xcode CLT, Windows needs Rtools.
- Install: `BiocManager::install("decontX", dependencies = TRUE)` or
  `devtools::install_dev_deps()` in a clone.
- First build compiles the Stan model and Rcpp code — expect several
  minutes; subsequent `devtools::load_all()` calls reuse the objects.

### Package-specific notes
- No Shiny app and no pkgdown site currently — the Shiny and pkgdown
  sections of the common playbook do not apply yet (see dev/ROADMAP.md).
- GENERATED files — never hand-edit: `R/RcppExports.R`,
  `src/RcppExports.cpp` (Rcpp), `R/stanmodels.R`, `src/stanExports_*`
  (rstantools). Edit `inst/stan/shrinkage.stan` and re-generate instead;
  changes to the Stan model are structural (ADR first).
- `decontPro()` runs MCMC sampling — tests and examples in that area are
  slow by nature; keep test fixtures tiny.
- Heavy Suggests (TENxPBMCData, SingleCellMultiModal, scran) are only
  needed for vignettes; `make check` runs with
  `_R_CHECK_FORCE_SUGGESTS_=false`.
- `R/celda_functions.R` mirrors helpers in the celda package — keep
  divergence minimal and coordinate cross-package changes via issues.
