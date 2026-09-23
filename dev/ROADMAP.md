# decontX Roadmap

Where the package is going. Proposals (from humans or agents) should
lean in these directions; anything structural still needs an ADR first.

Rewritten 2026-09-21 from the comprehensive review in
`dev/audits/2026-09-comprehensive-review.md` (finding IDs below refer to
it). Near-term work is organized as waves; each wave lands as its own
PR(s), with `make test` green throughout and NEWS.md entries for
user-facing changes.

## Near term

### Wave 1 — Make devel green + quick wins (no ADR needed)

- **celda dependency: keep it.** The deprecated-`normalizeCounts` warnings
  (D3-01, CI-01) that surface in decontX's examples originate in celda and
  are being fixed upstream in celda separately — decontX does not drop
  celda or work around the deprecation. No ADR-D1. (The dead `.cDCalcEM*`
  functions are still removed as ordinary dead-code cleanup below; that is
  independent of the dependency.)
- Fix the defunct Seurat call: `GetAssayData(slot=)` →
  `SeuratObject::LayerData(layer=)` (D3-03) + a mocked regression test.
- Fix the P0 numeric-`batch` indexing bug (`as.character(batch)`) with a
  regression test (D1-01).
- Delete dead code: commented-out legacy init block, 5 unreachable R
  internals, 6 unregistered C entry points + 3 dead Sparse variants,
  `nonzero`/`fastNormPropLog`/`fastNormPropSqrt` (keep `fastNormProp`),
  `R/package_skeleton_test.R`, permanently-unrunnable commented tests
  (D1-22..24, D2-04, D2-15, D4-16). This alone clears most of the
  commented_code lint backlog and ~300 uncovered lines.
- Metadata: URL/BugReports, expanded biocViews, inst/CITATION, NEWS
  backfill (drafts in `dev/audits/2026-09-review-drafts.md`); reconcile
  README.md/README.Rmd (D6-09); fix the stale `\link[uwot]` (D3-13).
- CI hygiene: checkout@v4 everywhere, pin upload-artifact@v4, add
  `_R_CHECK_FORCE_SUGGESTS_=false` to R-CMD-check, codecov.yml +
  `make coverage`, decide the oldrel-1 line, create `dev/agent-log.md`
  (D4-15, D5-02..06).
- Lint burn-down after the deletions (mostly indentation/return), then
  keep the lint job blocking.

### Wave 2 — Correctness hardening (no ADR)

- Move matrixSums registration out of generated `RcppExports.cpp` into a
  hand-owned `src/init.c` (D2-03) — **must precede anything that reruns
  `compileAttributes()`**.
- Input-validation pass: `.checkDecontXParams()` (maxIter/iterLogLik/
  convergence/seed), varGenes/dbscanEps check order, batch length,
  background gene-row matching (fixes issue #43's failure mode),
  bgBatch hoisting, NA/gap-level `z`, preserve user cluster labels,
  delta bounds, negative counts (D1-02..12).
- C-level guards: NA-group check in `_colSumByGroup_numeric` (+ its row
  twin if kept) — removes a reproducible segfault (D2-01); `min(z)`
  guard in `decontXInitialize` (D2-02); `R_xlen_t` offsets (D2-07);
  empty-cluster division guard in `decontXEM` (D2-08).
- `plotDensity` shared bandwidth (D7-04) and the decontaminated ≤
  observed invariant test (D7-05).

### Wave 3 — Test debt (guided by the 18-file plan in the drafts doc)

- Convert the two empty test blocks into asserted tests; assert the
  batch path; add the ground-truth recovery oracle (D4-01..03).
- Shrink the decontPro fixture (14.8s → 1.9s measured) and mock
  `.call_stan_vb` for Stan-free dispatch/validation tests (D4-05).
- Snapshot policy for user-facing messages; no-vdiffr plotting policy
  (structure/data assertions) — record both in AGENTS.md (D4-08, D4-14).
- Target: honest ~80% coverage with the suite wall clock halved.

### Wave 4 — Issue fixes and dependency diet (ADRs D2/D3/D5)

- **Issue #46:** default `output_samples = 10` in `.call_stan_vb` and
  expose `output_samples`/`seed`/`iter` on `decontPro()` (bit-identical
  results verified; small ADR since it adds args to an exported generic;
  merges D1-13/D6-08/D7-01). Document the TMPDIR footprint and the
  residual PSIS pass.
- **ADR-D2:** Seurat → SeuratObject in Imports (45 → 13 transitive).
- **ADR-D3:** drop plyr (`match(z, unique(z))`) and reshape2 (two small
  base-R helpers) — gated on Wave 3's plotting tests.
- **ADR-D5:** DESCRIPTION floors (drop Matrix pin, verify scrapper pin,
  R >= 4.5.0 on devel); drop dead scran Suggest; resolve the undeclared
  HDF5Array use (prefer deleting the redundant special case).
- Answer/close the open issues: #40 (corrected explanation + FAQ), #39
  (stale library), #43/#44/#45 (FAQ + the Wave-2 background validation).

### Wave 5 — Structure (ADR first)

- Consolidate the decontX argument surface (~11 duplicated sites with
  real default drift) into one params constructor + delegating methods
  (D1-20). Decompose `.decontX`/`.decontXoneBatch` into named helpers
  (D1-21) and fix the per-batch same-seed RNG reuse (D1-16).
- Migrate `Eigen::MappedSparseMatrix` → `Map<SparseMatrix<double>>`
  (D2-10/D3-15) — after the Wave-2 registration move.
- Regenerate the committed Stan exports with current rstantools (the
  committed files are stale — B2-01); note the GPL-3 template header in
  LICENSE/NOTICE (D2-13).
- Resolve the docs-vs-code question on always-computed UMAP (D6-04):
  either document reality or skip init when `z` is supplied.

### Quality baseline (carried over, now concrete)

- Get `make bioccheck` clean and flip the CI BiocCheck job to required
  (Wave 1 clears most NOTEs; function-length and indentation NOTEs
  shrink with Waves 1/5).
- Docs: FAQ sections in both vignettes (drafts ready — the #40 answer
  must use the corrected D7-03 version), `@return` expansions, decontPro
  reproducibility documentation, 1-based `cell_type` warning, internal
  man pages `@noRd`.

## Medium term

- **pkgdown site**: unchanged plan — add `_pkgdown.yml`, build locally,
  deploy to `gh-pages` per the lab playbook (CI verifies structure with
  `pkgdown::check_pkgdown()` only).
- **Precomputed vignettes**: adopt the `.Rmd.orig` pattern (decontPro
  certainly; decontX arguably — it currently runs the EM twice per
  build); interim: `eval=FALSE` the second decontX run (D6-11). Record
  as an ADR when done.
- **Community surface**: issue templates (sessionInfo + reprex), CODE_OF
  _CONDUCT link, NEWS checkbox in the PR template (D5-05, D6-14, D6-02).
- **Expand the decontPro vignette**: it is bare-bones today. Beyond the
  Wave-1 fixes (D6-13), grow it toward parity with the decontX vignette:
  what the three output matrices mean and how to access them, prior
  tuning guidance (`delta_sd`/`background_sd`) with visual before/after,
  the `ambient_counts` empty-droplet workflow, cluster-granularity
  advice (from the issue #40 discussion), reproducibility/`output_
  samples` notes once Wave 4 lands, and the corrected #40 FAQ.

## Longer term / ideas

- Performance work on the DecontX EM loop guided by profiling, not
  speculation — candidates already identified: the per-iteration C++→R
  `fit_dirichlet` callback (D3-07), the exp/log round-trip in
  `calculateNativeMatrix` (D2-09), and the per-batch matrix copies
  (D1-27).
- decontPro scalability beyond `output_samples`: the model has 2·N·M
  latent parameters, so memory stays linear in cells even with one draw;
  hierarchical/low-rank background structure is the real lever —
  methods-level, needs a methods review first (D7-02).
- **Native Seurat and AnnData support**: S4 methods/accessors so
  `decontX()` and `decontPro()` run directly on Seurat and AnnData-backed
  objects and write results back in place, sparing users the manual
  SCE round-trip. Today only `decontPro` has a Seurat method (input
  only); `decontX` has none, and neither supports AnnData (likely via
  zellkonverter/anndata — dependency choice is part of the design).
  New exported methods on the generics → ADR first; builds on the
  Wave-4 SeuratObject migration (ADR-D2).
- Keep `R/celda_functions.R` in sync with celda; coordinate upstream so
  celda re-exports decontX's functions instead of shipping a 5-year-old
  duplicate under the same names (D3-02) — ADR territory.

## Non-goals (for now)

- No Shiny app in this package (interactive use is served by
  singleCellTK).
- No change to the statistical models without a methods-level review —
  model changes are science, not refactoring. (This includes the
  decontPro low-rank idea above.)
