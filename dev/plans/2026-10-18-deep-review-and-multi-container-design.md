# Deep review, multi-container support, and docs: phased plan

- **Date:** 2026-10-18
- **Status:** Approved by the maintainer, 2026-10-18
- **Base:** `devel` at 3d7a8aa (decontX 1.11.2)
- **Builds on:** `dev/audits/2026-09-comprehensive-review.md` and the
  waves in `dev/ROADMAP.md`

This is the design for a multi-phase program. It is the spec the work is
built from; there is no separate step-by-step implementation plan. Each
phase still has its own gates (ADRs, approvals) listed in
[Sequencing and gates](#sequencing-and-gates).

## 1. Intent

### What the maintainer asked for

- **Phase I:** a deep review of the repository to find bugs, increase
  test coverage, improve code quality, and make the code faster and more
  memory-efficient; then fix what it finds.
- **Phase II:** make `decontX()` and `decontPro()` run easily on
  SingleCellExperiment, ExperimentSubset (the
  [ExperimentSubset](https://bioconductor.org/packages/ExperimentSubset)
  package's subset classes), Seurat, and AnnData objects.
- Review the vignettes and tutorials for simplicity, clarity, and
  completeness, especially decontPro.
- A phased plan for all of it.

### Decisions made during design

| Topic | Decision |
|---|---|
| Relation to the 2026-09 review | Re-baseline and delta-review; don't start over. |
| Review process | Sequential subagents, one area at a time, steerable. |
| "Easily" in Phase II | Same class in, same class out, results written back into the user's object, for every supported class. |
| AnnData user | R users holding `.h5ad` files. Python/scanpy users are out of scope. |
| AnnData package | anndataR preferred; the Phase II ADR confirms it is the better-maintained choice (vs zellkonverter). |
| Phasing | Review once up front; fixes; the core refactor joins Phase I and Phase II; docs in two passes. |
| Benchmarks | Simulated data for most work; a final benchmark on real PBMC data. Script committed under `dev/bench/`. |
| Coverage target | >= 85% on `R/`, >= 80% overall, every exported method and every `stop()` path tested. |
| Refactor and speed work | Bit-identical results, or a documented tolerance where floating-point reordering is unavoidable; anything beyond that goes back to the maintainer. |
| `decontPro()` return type | **Clean break**: container methods return the updated object. NEWS marks it Breaking. |
| Seurat split layers | Joined in memory; the layer split becomes `batch` when `batch` is `NULL`. |
| `z` / `batch` | Accept a vector or the name of a column in the object's cell metadata, for every class. |
| pkgdown | In scope, deployed to GitHub Pages. No development mode: the site documents the current release. |
| Vignette and article rebuilds | Local only, re-run at each release. No CI rendering. |
| Bioc vignettes | Short and fast; link to full-data articles once the site is live. |

### Success criteria

- One audit that every later phase works from, with every finding traced
  to a phase.
- Phase I: all confirmed P0/P1 bugs fixed with regression tests; coverage
  targets met; measured speed and memory improvements with before/after
  numbers; no unintended change in results.
- Phase II: each supported class runs `decontX()` and `decontPro()` in
  one call and gets its results back in the object, with parity tests
  proving the results equal the matrix method's.
- Phase III: short, fast Bioc vignettes a novice can follow as pasted; a
  complete decontPro vignette; a pkgdown site with full-data articles.

## 2. Current state (verified 2026-10-18)

- `.decontX()` (`R/decon.R:354`) and `.decontPro()` (`R/decontPro.R:125`)
  already take plain matrices. Container coupling sits in:
  - the exported methods (`R/decon.R:140-228`, `R/decontPro.R:64-122`);
  - `.decontxInitializeZLegacy()`, which builds an SCE internally (fine:
    internal, and needed by scater);
  - the plotting functions, which read `metadata(x)$decontX` and
    `colData(x)` directly (`R/plot_decontx.R:31`, `:467-477`).
- `decontPro()` returns a plain list for every input class. Its SCE method
  hardcodes the `"counts"` assay and has no `altExp` support.
- **Bug found during design:** on a Seurat object whose counts are split
  into layers (`counts.A`, `counts.B`), `LayerData(seu, layer = "counts")`
  returns only the first layer with a warning (verified with SeuratObject
  5.4.0). `decontPro()` on such an object silently drops every sample but
  the first.
- **Doc defects found during design:** the decontX vignette's Seurat
  section uses the defunct `GetAssayData(slot = )` and mixes
  `seuratObject`/`seuratObj`; the decontPro vignette uses `%>%` pipes and
  dplyr, explains only `decontaminated_counts`, and has no
  `ambient_counts` example or FAQ.
- 2026-09 waves still open in part: Wave 2 (input validation, C guards),
  Wave 4 (plyr/reshape2, DESCRIPTION floors), Wave 5 (params
  consolidation, `.decontX` decomposition, `MappedSparseMatrix`).
- ExperimentSubset classes inherit from their input class (a
  `SubsetSingleCellExperiment` is a SingleCellExperiment), so
  `decontX(es)` probably dispatches to the SCE method today but writes to
  the top level, not the subset. Not yet verified (ExperimentSubset isn't
  installed locally).
- anndataR's AnnData objects are R6 (`AbstractAnnData` base, with
  `InMemoryAnnData`, `HDF5AnnData`, `ReticulateAnnData`), cells x genes,
  with `X`, `layers`, `obs`, `var`, `obsm`, `uns`.

## 3. Phase 0: the review

Read-only: it changes no package code. The only commit it makes to the
package tree is the benchmark script under `dev/bench/`.

### Baseline

On a fresh `devel`: `make check`, `make bioccheck`, `make lint`,
`make test`, `make coverage` (per-file numbers). Record the status of
every 2026-09 finding: fixed, open, or obsolete.

### Review areas (one subagent each, in this order)

1. **R code: bugs.** Reproduce every claim live. Confirm the open
   2026-09 items (input validation, cluster-label preservation,
   background gene-row matching for #43, per-batch seed reuse) and the
   Seurat split-layer bug above. Check whether `z`/`batch` accept a
   column name today.
2. **Compiled code.** `src/`, including D2-01/02/07/08 and the
   `MappedSparseMatrix` deprecation.
3. **Tests and coverage.** Uncovered paths, hollow tests, missing error
   paths and edge cases; a per-file plan to reach the coverage targets.
4. **Code quality.** Duplication (the ~11 copies of the argument
   surface), function size, naming, the plyr/reshape2 removal.
5. **Speed and memory.** Profile `decontX()` with profvis and
   `bench::mark()` on simulated data at about 10k, 50k and 100k cells,
   plus one batched run; profile `decontPro()` on a realistic ADT panel.
   Rank hotspots by measured cost. Known candidates: the `fit_dirichlet`
   callback, the exp/log round-trip in `calculateNativeMatrix`, per-batch
   matrix copies.
6. **Object coupling.** Every place code touches a container directly;
   what each target class needs for read and write-back; whether
   ExperimentSubset objects dispatch to the SCE method today and what
   goes wrong.
7. **Docs.** Both vignettes, the README, the shared example
   (`man/examples/decontX.R`): clarity, completeness, the novice style
   rules, and correctness. Measure each vignette's build time and
   downloads. Dedicated decontPro gap analysis. Settle the DecontPro
   citation year (README says 2024, AGENTS.md says 2023).
8. **Targeted methods check.** Open usage issues #40, #42, #43, #44, #45,
   only where they point at code or docs.

### Benchmark script: `dev/bench/`

- `dev/bench/bench-decontx.R`: simulated data from
  `simulateContamination()` at fixed seeds and sizes; reports wall time
  and peak memory for `decontX()` (single and batched) and `decontPro()`.
  A `real` mode runs the final PBMC benchmark (TENxPBMCData).
- Uses bench and profvis, which stay dev-only tools, not DESCRIPTION
  dependencies. `dev/` is already in `.Rbuildignore`.
- Outputs go to the scratchpad; numbers are copied into the audit and PR
  descriptions.

### Deliverables

- `dev/audits/2026-10-deep-review.md`: findings with ID, severity
  (P0-P3 as in 2026-09), effort, evidence, and target phase.
- `dev/ROADMAP.md` rewritten to this plan's phases.
- A list of proposed GitHub issues and ADR stubs. Nothing is filed until
  the maintainer approves.
- An entry in `dev/agent-log.md`.

Gotchas: covr and bioccheck regenerate `src/stanExports_shrinkage.*`
(restore with `git restore src/`) and leave `decontX.BiocCheck/` and a
tarball in the repo root.

## 4. Phase I: fixes, coverage, quality, speed

A series of small PRs into `devel`, one issue per PR, each approved
before the next starts. Every PR bumps z, adds NEWS for user-facing
changes, keeps `make test` green, and never lowers coverage.

### I-a: regression tests

Before a fix lands, add tests for the behavior Phase 0 marks untested in
the code it touches, using `simulate_oracle()` and snapshots for
messages. These tests pin current behavior so later changes are
measurable.

### I-b: correctness fixes

P0, then P1, then P2. Each starts with a failing regression test. A fix
that changes numeric results (for example the per-batch seed reuse) says
so in the PR and in NEWS. The Seurat split-layer bug is fixed here for
the existing `decontPro()` Seurat method; Phase II rebuilds that method
anyway, so the fix is the minimal one (join layers before reading).

### I-c: coverage

Work through the Phase 0 per-file plan to >= 85% on `R/` and >= 80%
overall, with every exported method and `stop()` path tested. Hollow
tests are tracked and fixed alongside the percentage. No fast/slow test
split.

### I-d: core refactor (ADR-0004)

ADR-0004, "Container-agnostic core", is approved before code. It defines:

- **The core contract.** `.decontX()` takes counts (genes x cells,
  matrix or sparse), `z`, `batch`, background counts, background batch,
  and a params object; it returns a plain list (`decontXcounts`,
  `contamination`, `z`, `estimates`, `runParams`, UMAP). `.decontPro()`
  likewise. The core never reads or writes a container.
- **One params constructor/validator** (`.decontXParams()`) replacing
  the duplicated argument surface and its default drift
  (D1-20). Method signatures keep their user-facing defaults.
- **Decomposition** of `.decontX()` / `.decontXoneBatch()` into named
  helpers in the same file (D1-21).

Bit-identical: tests compare old and new output on fixed seeds across
the single-batch, batched, background, and user-`z` paths.

Separate PRs in this sub-phase: drop plyr and reshape2; move matrixSums
registration to `src/init.c` and migrate `MappedSparseMatrix`; the lint
backlog (its own PR, per the standards).

### I-e: speed and memory

One optimization per PR, each with before/after numbers from
`dev/bench/` in the PR description. Results bit-identical, or within a
documented tolerance if floating-point reordering is unavoidable; larger
differences go back to the maintainer. The real-PBMC benchmark runs once
at the end of I-e and is recorded in the audit.

No changes to the statistical models (ROADMAP non-goal).

## 5. Phase II: multi-container support

### Architecture

Three layers:

1. **Core** (from I-d): matrix in, plain list out.
2. **One method per class** on `decontX()`, `decontPro()`, and
   `decontXcounts()` / `decontXcounts<-`. Each reads inputs out of the
   object, calls the core, writes results back. Methods stay thin; logic
   in a method means the core is missing something.
3. **An internal results getter per class**, used by
   `plotDecontXContamination()`, `plotDecontXMarkerPercentage()`, and
   `plotDecontXMarkerExpression()` so they accept any supported object.

### Read and write locations

| Class | Input | Decontaminated counts | Per-cell results | UMAP | Run info |
|---|---|---|---|---|---|
| SCE (SpatialExperiment by inheritance) | `assay(x, assayName)`; decontPro adds `altExpName` | assay `decontXcounts` (unchanged) | `colData` | `reducedDim(x, "decontX_UMAP")` | `metadata(x)$decontX` |
| ExperimentSubset | `assay(es, subsetName)` | `setSubsetAssay()` | `subsetColData<-` | `reducedDim(es, "decontX_UMAP", subsetName =)` | confirmed in the spike |
| Seurat | `LayerData(seu, assay =, layer = "counts")`, split layers joined | new assay `decontXcounts` (one joined layer) | `meta.data` | reduction `decontX_UMAP` (`CreateDimReducObject`) | `Misc(seu, "decontX")` |
| AnnData (anndataR) | `X` or `layers[[layer]]`, transposed to genes x cells | `layers[["decontXcounts"]]` | `obs` | `obsm[["X_decontX_umap"]]` | `uns$decontX` |

decontPro writes its decontaminated counts to the same kind of location
(assay or layer named `decontPro`, final name set in ADR-0005) and its
other outputs (ambient and background matrices, parameters) to the
run-info location.

### Conventions across classes

- `z` and `batch` accept a vector or a single string naming a column in
  `colData` / `meta.data` / `obs` (or the subset's colData).
- Accessors return genes x cells for every class.
- **Seurat split layers:** the method joins the count layers in memory,
  in `Cells(seu)` order. A supplied `batch` is used as-is; if `batch` is
  `NULL` and the layers are split, the layer split becomes the batch and
  a message says so.
- **AnnData is R6, so it's modified in place.** Docs show
  `adata <- decontX(adata)` and warn that the original object changes
  too. Writes to a backed `HDF5AnnData` go to disk; the method warns.
- **decontPro clean break:** the SCE, Seurat, ExperimentSubset, and
  AnnData methods return the updated object. The `ANY` (matrix) method
  keeps returning the list.

### Dependencies and dispatch (decided in ADR-0005)

- anndataR and ExperimentSubset go in **Suggests**, used behind
  `requireNamespace()` guards with an error naming the missing package.
  SeuratObject is already in Imports.
- AnnData methods dispatch on `setOldClass("AbstractAnnData")`, which
  needs no import.
- ExperimentSubset: its S4 classes can't take methods from a Suggests
  package cleanly, so the SCE method gains `subsetName = NULL` and
  branches on `methods::is(x, "SubsetSingleCellExperiment")`. ADR-0005
  records the alternative (ExperimentSubset in Imports, which pulls in
  TreeSummarizedExperiment) and why it was rejected or chosen.
- ADR-0005 also confirms anndataR over zellkonverter (recent commits,
  issue response, Bioc check status).

### Spike before ADR-0005

Throwaway code in the scratchpad, results reported to the maintainer:

1. S4 dispatch on an R6 subclass (`InMemoryAnnData`) via
   `setOldClass("AbstractAnnData")`.
2. How ExperimentSubset names subset assays, whether names must be
   unique across subsets, and whether subsets hold metadata.
3. Seurat split-layer joining and write-back on SeuratObject 5.x.
4. Whether nested `uns$decontX` (per-batch estimates) survives
   `write_h5ad()` / `read_h5ad()`.

### Testing

- **Parity tests** per class: same data and seed give results identical
  to the matrix method, for `decontX()` and `decontPro()`.
- Write-back tests: every location in the table is populated and
  readable through the accessors.
- An h5ad write/read round trip.
- Seurat split-layer tests (with and without `batch`).
- Error paths: missing optional package, unknown column name, unknown
  `subsetName`, wrong orientation.
- `skip_if_not_installed()` for ExperimentSubset, anndataR, rhdf5;
  decontPro fixtures stay at 10 x 8.

## 6. Phase III: docs and site

All user-facing code follows the novice style rules: one operation per
line, defaults only, no `pkg::` in code shown, `set.seed()` before
anything random, ending with `sessionInfo()`. Every edited vignette or
article is rendered locally with `make article FILTER=<name>` and gets a
**novice-reader pass**: a fresh subagent with no context runs it as
written and reports confusing steps and unexplained terms before the
maintainer reviews.

### III-a: clarity and completeness (parallel with Phase I)

- **decontPro vignette, rewritten**, short and fast:
  1. What DecontPro models (native, ambient, background), one paragraph
     and a small diagram.
  2. Loading data; choosing `cell_type` and cluster granularity (from
     the #40 discussion); clustering one operation per line, no pipes.
  3. Running it; what every output element means and how to access it.
  4. Tuning the priors: before/after for two `delta_sd` /
     `background_sd` settings, read out with `plotDensity()`.
  5. The `ambient_counts` empty-droplet workflow (bundled or simulated
     empty droplets if cord_blood doesn't provide them).
  6. Runtime, memory (#46), and reproducibility.
  7. FAQ.
- **decontX vignette:** fix the Seurat section (temporary until III-b);
  add an FAQ (the corrected #40 answer, #43 background, #44 doublets,
  #45 many samples via `batch`); stop running the EM twice where it
  isn't needed.
- **README:** sync with the vignettes; fix the citation year.
- **Data and build time:** Phase 0's measurements decide whether the
  decontX vignette keeps a subsampled PBMC4k or moves to simulated data,
  and how far the decontPro subsample shrinks. Short vignettes still run
  real code (BiocCheck flags heavy `eval = FALSE` use).

### III-b: multi-container docs (after Phase II)

- One short "Using your object" section per class in each vignette
  (SCE, ExperimentSubset, Seurat, AnnData/h5ad), replacing the manual
  round trip. decontPro covers `altExpName` for SCE and the ADT assay
  for Seurat.
- Runnable examples for each new method, guarded with
  `requireNamespace()` where the package is optional; `\donttest{}`
  where slow.
- NEWS entry and a vignette callout for the `decontPro()` break.

### III-c: pkgdown site and articles (ADR first)

The pkgdown ADR records the layout, the release-only rebuild policy, and
the deploy trigger.

- `_pkgdown.yml` with a reference index; **no development mode**, so the
  site always documents the current Bioc release.
- Full-data articles in `vignettes/articles/`: full PBMC for decontX,
  full cord_blood CITE-seq for decontPro. Heavy articles may use
  precomputed outputs.
- Articles are run and checked locally only. `dev/RELEASE.md` gains a
  step to re-run the articles and rebuild the site at each release.
- A gh-pages deploy workflow from the shared standards, triggered at
  release. Deploys go to the `gh-pages` branch; the `docs/` deny rule is
  unaffected. CI keeps only `pkgdown::check_pkgdown()`, and
  `make site-check` runs before PRs that touch the site.
- Site URL added to DESCRIPTION `URL:` and the README.
- Enabling GitHub Pages and the first deploy are maintainer actions.
- **After the first deploy**, the short vignettes link to the matching
  articles. The FAQ stays in the vignettes (it ships with the package);
  articles link back to it rather than repeating it.

## 7. Sequencing and gates

```
Phase 0 review ─► I-a tests ─► I-b fixes ─► I-c coverage ─► I-d core (ADR-0004) ─► I-e speed
     │                                                           │
     └─► III-a docs (parallel)            Phase II spike ─► ADR-0005 ─► II methods ─► III-b docs
                                                                                       │
                                     pkgdown ADR ─► III-c site + articles ─► first deploy ─► vignette links
```

Gates, each needing explicit maintainer approval:

1. Phase 0 audit and ROADMAP rewrite, before any issue is filed.
2. ADR-0004, before I-d.
3. Phase II spike results, then ADR-0005, before Phase II code.
4. The pkgdown ADR, before III-c.
5. Every PR, reviewed locally before it is pushed.

Release timing: the decontPro break and the new container methods land
in the same devel cycle, so one Bioc release carries the whole API
change. Phase I fixes can ship in any cycle; any P0 is flagged for a
possible release-branch port, which the maintainer decides.

## 8. Effects on related packages

- **celda** forwards to `decontX::decontX()` on SCE. That path stays
  bit-identical through I-d and Phase II, enforced by the parity and
  refactor tests. The column-name convention for `z`/`batch` only adds
  behavior.
- **singleCellTK** reaches decontX through celda; same as above.
- **ExperimentSubset** (lab package): gaps found in the spike (subset
  metadata, naming) are raised as issues there, not worked around here.
- **ADR-0003 Stage 2** (inverting the celda dependency) stays separate.

## 9. Out of scope

- Changes to the DecontX or DecontPro statistical models, including the
  decontPro low-rank idea.
- Python/scanpy users; MuData.
- A `SummarizedExperiment` method beyond what ExperimentSubset needs.
- CI rendering of vignettes or articles.
- A pkgdown development (`/dev/`) site.

## 10. Amendments (2026-10-27, from the Phase 0 review)

Decided by the maintainer during the review; details and evidence in
`dev/audits/2026-10-deep-review.md`. Where they conflict with the
sections above, these win.

- **ExperimentSubset (section 5):** the method uses the direct
  approach. It reads `assay(x, subsetName)` and the subset's colData,
  calls the core, and writes back with `setSubsetAssay()`,
  `colData(x, subsetName =) <-`, `reducedDim(x, , subsetName =) <-` and
  `metadata(x, subsetName =) <-` (A6-12). It does not use
  `getSubsetExperiment()`/`setSubsetExperiment()`. decontX adds no
  ExperimentSubset detection code (no `is(x,
  "SubsetSingleCellExperiment")` branch, which also misses two of the
  three classes) until ExperimentSubset #5 lands. The Suggests floor is
  the ExperimentSubset release with #53 (#54 isn't needed). Gaps found
  were filed upstream: ExperimentSubset #2-#5 and #57-#59.
- **AnnData (section 5):** `setOldClass("AbstractAnnData")` doesn't
  dispatch (A6-03). The ANY methods redirect on
  `inherits(x, "AbstractAnnData")`, and the AnnData method extracts the
  matrix and calls the core directly, with no temporary SCE.
- **decontPro (sections 5, 6):** the empty-droplet argument is renamed
  `background`, matching decontX, with `ambient_counts` kept as a
  permanent alias. The `ambient_counts` empty-droplet workflow (III-a,
  item 5) is dropped from the vignette and the full-data article;
  `background` gets a short example like decontX's. decontPro's argument
  changes get their own ADR (ADR-0006).
- **HDF5 input (Phase I-d):** returns an in-memory sparse result wrapped
  as a DelayedMatrix, instead of writing back to a temporary HDF5 file
  (A5-06). NEWS.
- **Initialization (Phase I-e):** runs by default when `z` is given, as
  today, with a new opt-out argument; plus opt-in `numThreads` (A5-01,
  A5-02).
- **ADR numbering:** ADR-0004 core, ADR-0005 container methods,
  ADR-0006 decontPro arguments, ADR-0007 pkgdown.
