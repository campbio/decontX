# decontX Roadmap

Where the package is going. Proposals (from humans or agents) should
lean in these directions; anything structural still needs an ADR first.

Rewritten 2026-10-27 from the deep review in
`dev/audits/2026-10-deep-review.md` (finding IDs below refer to it), and
organized by the phases of
`dev/plans/2026-10-18-deep-review-and-multi-container-design.md`. It
replaces the 2026-09 waves; their open items are folded in. Every PR
keeps `make test` green, never lowers coverage, bumps z, and adds NEWS
for user-facing changes. Each PR is reviewed locally before it is
pushed, one at a time.

## Done

- **Phase 0 (review), approved 2026-10-27:** baseline, status of the
  2026-09 findings, eight review areas, the benchmark script
  `dev/bench/bench-decontx.R`, and this rewrite. Upstream
  ExperimentSubset issues filed: #2-#5, #57-#59. The work below was
  opened as #65-#84 (2026-10-29); replies posted on #40, #43, #44, #46,
  #39 and #56, ADR-0004 to ADR-0007 drafted as Proposed.
- **Phase I-a to I-c (1.11.3 to 1.11.13), 2026-10-30 to 2026-11-04,
  on a devel (Bioc 3.24) library set through an untracked `.Renviron`:**
  - Stopgap: the decontX vignette's Seurat section (#65, PR #86; A7-02).
  - I-a: stored reference outputs and `dev/compare-commits.R` (#66,
    PR #87; A3-01 and the pins).
  - I-b: input validation (#67, PR #88); `background` and
    `ambient_counts` matched by row name, `bgBatch` checked (#68,
    PR #89); one cluster with background, k-means on PCs (#69, PR #90);
    user cluster labels and marker plots (#70, PR #91); small batches
    and input classes (#71, PR #92); decontPro inputs (#72, PR #93);
    plots (#73, PR #94); cosmetic (#74, PR #95). Follow-up: #96.
  - I-c: coverage 92.5% overall and 98.4% on `R/`, every exported
    method and `stop()` path tested, a warning-free and quiet suite
    (#75, PR #97). The `expect_equivalent()` warning (#55, A3-07) was
    fixed by an outside contribution, #64.
  - #43 closed once #68 and #69 merged.
- From the 2026-09 waves: green devel CI, celda removed from Imports
  (ADR-0003 Stage 1), the numeric-`batch` fix, the defunct Seurat call,
  the #46 memory fix (1.11.1), metadata and CITATION, lint backlog
  cleared.

## Phase I: fixes, coverage, quality, speed

### I-d: core refactor (ADR-0004 first), #76-#78

- **ADR-0004, container-agnostic core:** `.decontX()` and `.decontPro()`
  take matrices and return plain lists (decontaminated counts,
  contamination, labels, batch, estimates, UMAP, `runParams` with `seed`
  and `iterLogLik`); one `.decontXParams()` constructor with no
  defaults, full method signatures, and a test that the methods' formals
  agree; `.decontX` and `.decontXoneBatch` split into named helpers.
  Bit-identical on the reference paths (A4-01 to A4-08, A7-22).
- Batched assembly with one `cbind` and one reorder (A5-05); HDF5 input
  returns an in-memory DelayedMatrix (A5-06, NEWS).
- #76 lands as five PRs from and into devel (decided 2026-11-03; ADR-0004
  "Decided 2026-11-03"): (1) the core contract, `.decontXParams()`, the
  formals test, and #96; (2) decomposition; (3) batched assembly;
  (4) HDF5; (5) the file move. All five are done (PRs #100, #101, #103,
  #104 and the file move; 1.11.14 to 1.11.18).
- Separate PRs:
  - plyr replaced exactly, reshape2 by a faithful helper (#77; A4-09,
    A4-10). Done in 1.11.24, with ADR-0008;
  - dead compiled code deleted (#78), `.colSumByGroup` becomes `rowsum()`,
    `matrixSums.c` and `R/matrixSums.R` go (coordinate with celda), then
    the `MappedSparseMatrix` migration (D2-01, D2-04 to D2-07, D2-10,
    D2-11, D2-13 to D2-15, A2-02, A7-16). Done in 1.11.25 (spec
    `dev/plans/2026-11-09-dead-compiled-code.md`); celda keeps its own
    copies of the matrix sums, so nothing there changes;
  - scrapper >= 1.5 and the old-API branch removed; unreachable code
    removed; explicit `decontXcounts` methods (A3-11 to A3-13);
  - small quality items (A4-11, A4-13 to A4-16).

### I-e: speed and memory

One optimization per PR, with before/after numbers from `dev/bench/` in
the PR (real PBMC mode for EM changes). Bit-identical, or the documented
tolerance.
- Opt-out of initialization when `z` is given, and opt-in `numThreads`
  (#79; A5-01, A5-02). New arguments: in ADR-0004 or a short ADR.
  Done in 1.11.19 (`embedding`, `numThreads`, and `reducedDimName` in
  `plotDecontXContamination()`; spec
  `dev/plans/2026-11-07-embedding-numthreads.md`).
- Unchecked indexing in the EM loop (#80; A5-03); the direct
  `calculateNativeMatrix` (A5-04, D2-09, within 1e-14);
  `.process_stan_vb_out` tidy-up (A4-20). Unchecked indexing done in
  1.11.20 (spec `dev/plans/2026-11-08-em-unchecked-indexing.md`); the
  direct native counts in 1.11.21 (`calculateNativeCounts()`; spec
  `dev/plans/2026-11-08-native-matrix-direct.md`); the tidy-up in 1.11.22
  (spec `dev/plans/2026-11-08-stan-process-tidy.md`).
- The real-PBMC benchmark once at the end, recorded in the audit.

## Phase II: multi-container support (spike results, then ADR-0005)

decontPro's arguments change in the same devel cycle, so one Bioc release
carries the whole API change.
- **ADR-0006, decontPro arguments (#81):** camelCase primary with permanent
  snake_case aliases (error if both given); the empty-droplet argument
  becomes `background` with `ambient_counts` as an alias; `iter`, `seed`,
  `output_samples`; a quiet option; the generic dispatches on the first
  argument only; rstan's temporary CSV removed when the fit ends (A5-08,
  A4-12, A4-14).
- **ADR-0005, container methods (#82):** same class in and out, results
  written back, `z`/`batch` as a vector or a cell-metadata column name,
  and decontPro container methods return the object (clean break, NEWS
  marks it Breaking). Write locations are the audit's area 6 table.
  - Seurat: join split layers on a copy; the layer split becomes `batch`
    when `batch` is NULL; separate `decontXcounts` assay; per-batch
    reductions hold only that batch's cells (A1-05, A6-08).
  - AnnData via anndataR (Suggests): extract the matrix, call the core,
    write to `layers`/`obs`/`obsm`/`uns`; dispatch by an `inherits()`
    redirect in the ANY methods; clear errors for read-only backed
    objects (A6-03, A6-10). Restore the names and dimnames `uns`
    drops, including UMAP column names, when rebuilding the estimates
    (A6-04; #73's positional plot fix is a stopgap).
  - SCE: decontPro reads and writes an altExp (`altExpName`).
  - Plots accept any supported object through an internal getter
    (A4-18, A4-19).
  - ExperimentSubset last, after ExperimentSubset #53 and #5 are in a
    release: the direct approach (A6-12 to A6-19), Suggests floor on
    that release, no workaround code.
- Parity tests against the matrix method for every class; write-back
  tests; an h5ad round trip; Seurat split-layer tests.
- SeuratObject to Suggests is a low-priority option (A5-09).

## Phase III: docs and site

### III-a: clarity and completeness (alongside Phase I), #83

- **decontPro vignette rewrite** (A7-03 to A7-05, A7-11): the model;
  cord_blood loaded with only the needed modes; its bundled cell-type
  labels (no Seurat, no dplyr, which leaves Suggests); 500 cells;
  what every output means; tuning the priors with before/after plots; a
  short `background` example like decontX's; runtime and
  reproducibility; FAQ. No empty-droplet demo.
- **decontX vignette:** show the cell-type mapping code (A7-01);
  "proportion", not "percentage" (A7-08); simpler code (A7-09); one EM
  run (A7-06); FAQ with #40, #43, #44 (decontX and doublet detectors on
  the same original counts), #45 (A7-07, A8-03).
- Citations in both vignettes and the package help page; 2024 for
  DecontPro everywhere, including `AGENTS.md` (A7-12, A7-14).
- Roxygen and examples (A7-10, A7-15, A7-18, A7-23); README (A7-20);
  DESCRIPTION Title and Description (A7-21); chunk labels, long lines,
  typos (A7-19, #57).
- Every edited vignette is rendered locally and gets a novice-reader
  pass.

### III-b: multi-container docs (after Phase II)

- A short "Using your object" section per class in each vignette, which
  replaces the manual round trip; runnable examples for each method;
  NEWS and a vignette callout for the decontPro break.

### III-c: pkgdown site and articles (ADR-0007 first), #84

- Documentation moves to GitHub Pages: a release site at the root,
  rebuilt at each release from the `RELEASE_*` branch, and a development
  site at `/dev/`, built from `devel`. The pkgdown mode is set by branch
  (Bioc's three-part versions defeat `mode: auto`).
- `_pkgdown.yml`; full-data articles for PBMC (decontX) and cord_blood
  (decontPro), run locally at each release; how the articles appear on
  the development site is decided in ADR-0007. Enabling Pages and the
  first deploy are maintainer actions. After it, the vignettes link to
  the articles.

## Longer term

- **A spatial model** (Xenium and others, #42): a new method, planned by
  the maintainer after this program. Not part of it.
- **ADR-0003 Stage 2:** invert the celda/decontX relationship so celda
  uses decontX's functions, sequenced across both packages' Bioconductor
  releases so neither breaks for users. See
  `dev/adr/0003-remove-celda-dependency.md`.
- decontPro scalability beyond `output_samples`: memory stays linear in
  cells; hierarchical or low-rank background structure needs a methods
  review first.

## Non-goals (for now)

- No Shiny app in this package (interactive use is served by
  singleCellTK).
- No change to the DecontX or DecontPro statistical models, or to
  decontPro's VB defaults, without a methods-level review.
- Python/scanpy users and MuData.
- CI rendering of the full-data articles. (The development site renders
  the short vignettes in CI.)
