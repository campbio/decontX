# Deep review: decontX (Phase 0)

- **Date:** 2026-10-18 to 2026-10-27
- **Version reviewed:** 1.11.2, `devel` at c6fca8e (the area 6
  ExperimentSubset re-check ran at 867e2ca, which changes no package code)
- **Plan:** `dev/plans/2026-10-18-deep-review-and-multi-container-design.md`
  (Phase 0). This audit is the evidence; `dev/ROADMAP.md` is the
  actionable view and was rewritten from it.
- **Evidence:** `dev/audits/2026-10-review-drafts.md` holds each area's
  full notes: reproductions, measurements, line references and the
  decisions as they were made. The scripts they mention lived in a
  session scratchpad and are gone; the findings, inputs and outputs are
  written out in the drafts.
- **Method:** a baseline on `devel`, a status pass over the 2026-09
  findings, then eight review areas run as sequential subagents. Every
  finding was reproduced live unless marked "code read". The lead
  re-verified the claims marked in each area's heading. The audit changes
  no package code; the only package-tree commit was the benchmark script
  `dev/bench/bench-decontx.R` (10184ce).

Severity: P0 = wrong results or breakage users hit, P1 = real risk of
wrong results or a crash, P2 = worthwhile, P3 = nice to have. Effort:
S < half a day, M = 1-3 days, L = larger. Phase is the plan's phase:
I-a regression tests, I-b correctness fixes, I-c coverage, I-d core
refactor (ADR-0004), I-e speed and memory, II multi-container
(ADR-0005), III-a docs, III-b multi-container docs, III-c pkgdown site.

## Summary

1. **No P0 findings.** Devel CI is green, `make check` is clean, lint
   is clean, and 23 of the 2026-09 findings are fixed, including all of
   the P0s.
2. **Eleven P1 bugs, all silent or crashing input-handling failures.**
   The worst give wrong results without an error:
   - a `z` shorter than the number of cells runs, and pads with NA
     (A1-01);
   - `background` gene rows are never matched to the counts, so reordered
     genes give wrong estimates and extra genes give issue #43's error
     (A1-02); decontPro's `ambient_counts` has the same bug (A1-12);
   - a `bgBatch` that doesn't match drops the background silently
     (A1-03);
   - decontPro on a Seurat object with split layers runs on the first
     layer only, and decontPro ignores the ADT altExp (A1-04, A1-05);
   - the marker-percentage plot is wrong for cluster labels that aren't
     1..K (A1-07).
3. **The tests can't yet catch a change in results.** There are no
   stored reference outputs (A3-01), coverage of `R/` is 75% against an
   85% target, and 26 of 35 `stop()` paths are untested (A3-02).
   `plot_decontPro.R` has no coverage at all.
4. **Speed: about 90% of decontX's runtime is the initialization (PCA,
   kNN, UMAP), and it runs even when the user supplies `z`.** The EM is
   about 5% on simulated data and 33% on PBMC4k. The most effective
   change is an option to skip initialization when `z` is given (A5-01);
   the EM itself has one bit-identical win (A5-03). Batched runs make
   about 7.5 full copies of the counts (A5-05). HDF5 input spends 63% of
   its time writing the result back to a temporary HDF5 file (A5-06).
5. **decontPro is 98% Stan gradient time.** There are no bit-identical
   speedups; the model and VB defaults are out of scope.
6. **Multi-container support is feasible for every target class without
   temporary SCE objects.** Seurat (SeuratObject 5.4) and AnnData
   (anndataR 1.2.2) prototypes worked. ExperimentSubset works with the
   direct approach on the ExperimentSubset PR stack (#53); upstream fixes
   were filed instead of workarounds (ExperimentSubset #2-#5, #57-#59).
7. **The shipped vignettes break when pasted.** The decontX vignette
   hides the code that builds `cellTypeMappings` (A7-01), and its Seurat
   section fails three ways (A7-02). The decontPro vignette prints about
   45 lines of Stan warnings and explains one of its four outputs (A7-03).
8. **The 2026-09 doublet advice was wrong** (A8-03). The FAQ now says to
   run decontX and the doublet detector on the same original counts.

## Baseline

At c6fca8e, 2026-10-18, R 4.6.1, local Bioc 3.23 (release) with
BiocCheck 1.48.1. The package targets devel (3.24); Phase I sets up a
separate devel library first.

| Check | Result |
|---|---|
| `make lint` | no lints |
| `make test` | 84 pass, 0 fail, 2 warnings (`expect_equivalent`, #55; `CreateSeuratObject` coercion), 0 skips, about 33 s |
| `make check` | OK |
| `make bioccheck` | 1 ERROR (support site HTTP 504, transient), 1 WARNING (release-version check, environment), NOTES below |
| Coverage (covr) | 74.21% overall; 71.3% on codecov's denominator, which excludes generated Stan code |
| Devel CI | green since #63 (2026-10-09) |

Coverage by file: `plot_decontPro.R` 0%, `matrixSumsSparse.cpp` 0%,
`matrixNorm.cpp` 17%, `celda_functions.R` 30%, `matrixSums.R` 64%,
`decon.R` 84%, `plot_decontx.R` 85%, `DecontX.cpp` 86%, `decontPro.R`
96%, `matrixSums.c` 98%, `stan_helpers.R` 100%.

BiocCheck NOTEs that are defects: R dependency floor 4.3.0 (#57); missing
chunk labels in both vignettes (#57); `cat()` in
`R/celda_functions.R:14` (A4-08); 14 functions over 50 lines (A4-03);
`decontXcounts.Rd` has no examples (A7-15); long lines in vignettes
(A7-19). False positives: `sapply` in generated `R/stanmodels.R`, lines
not indented by 4 spaces (the package uses 2, per `.lintr`), generated
lines over 80 characters.

## Status of the 2026-09 findings

Of 75 findings in `dev/audits/2026-09-comprehensive-review.md`: 23 fixed,
11 partly fixed, 40 open, 1 not applicable (D5-04). The fixed ones
include every P0 (D1-01 numeric batch, D3-03 defunct Seurat call, celda
dependency) and the #46 memory fix. The open ones are folded into this
audit's findings, which reference the old ID where they overlap (for
example A1-02 is D1-03). The full list is in the drafts.

## Findings

### Area 1: R code bugs

| ID | 2026-09 | Finding | Sev | Eff | Phase |
|---|---|---|---|---|---|
| A1-01 | new | `z` (or `batch`) shorter than the number of cells runs silently; padding NA becomes its own cluster (correlation with truth -0.02 vs 0.999). `decon.R:499-501` | P1 | S | I-b |
| A1-02 | D1-03, #43 | `background` rows never matched to the counts: reversed rows change contamination by up to 0.47; extra rows give #43's error; non-overlapping names run. `decon.R:1268-1315` | P1 | S | I-b |
| A1-03 | D1-08 | `bgBatch` mismatch silently drops the background (correlation 0.988 -> 0.746); short vectors recycle. `decon.R:480, 1298-1306` | P1 | S | I-b |
| A1-04 | new | decontPro on Seurat v5 split layers reads only the first layer, with a warning. `decontPro.R:93, 95` | P1 | S | I-b (full fix II) |
| A1-05 | new | decontPro reads the default Seurat assay or SCE `counts`, so it ran on RNA, not ADT; the altExp is ignored. `decontPro.R:70-95` | P1 | M | II |
| A1-06 | new | A DelayedMatrix with a dense in-memory seed crashes after the full EM ("the condition has length > 1"). `decon.R:617-618` | P1 | S | I-b |
| A1-07 | new | `plotDecontXMarkerPercentage` is wrong for numeric labels that aren't 1..K (0-based errors, 10/11/12 give NA, 2/5/9 mislabelled; unused factor levels give Inf). `plot_decontx.R:225-229, 572-581` | P1 | S | I-b |
| A1-08 | D1-02 | A one-cell batch, or one-cell complement, crashes ("subscript out of bounds"). `decon.R:479, 480, 569` | P1 | S | I-b |
| A1-09 | D1-06 | User `z` labels come back as integer codes; factor codes single-batch but labels batched; NA crashes. `decon.R:453, 589-594` | P1 | S | I-b |
| A1-10 | D1-04/05 | Scalar parameters unvalidated: `maxIter = 0` cryptic, `convergence = "x"` stops after one iteration, `estimateDelta = "yes"` becomes FALSE, and more. `decon.R:779-822` | P1 | S | I-b, consolidated I-d |
| A1-11 | D1-13 | decontPro `cell_type`/counts unvalidated: gaps, NA, length; a Stan failure surfaces as "no applicable method for '@'". | P2 | S | I-b |
| A1-12 | new | decontPro `ambient_counts` rows not matched (decontPro's #43): reversed rows accepted, ambient profile wrong (max diff 0.93 on a real fit). `decontPro.R:137-147` | P1 | S | I-b |
| A1-13 | new | decontPro's `...` swallows unknown arguments (`seed`, a `delta_SD` typo, `batch`). | P2 | S | I-b |
| A1-14 | new | Dense input returns a dgCMatrix; the coercion back is dead code. | P2 | S | I-b |
| A1-15 | new | dgRMatrix errors, dgTMatrix warns, data.frame input is cryptic. `decon.R:491, 495` | P2 | S | I-b |
| A1-16 | D1-15 | Initialization runs when `z` is supplied; the zero-count message tells users to supply `z` when they already did. | P2 | S | I-b (message), I-e (A5-01) |
| A1-17 | D1-11 | `delta = c(0, 0)` accepted; NA delta cryptic. | P2 | S | I-b |
| A1-18 | D1-12 | Negative counts pass through to the output. | P2 | S | I-b |
| A1-19 | D1-09 | `batch` length and NA unvalidated. | P2 | S | I-b |
| A1-20 | new | Marker plots with one marker, or permuted numeric-like colnames, fail or give NA clusters. `plot_decontx.R:214, 391-400` | P2 | S | I-b |
| A1-21 | new | `simulateContamination()` ignores `NRange` when min equals max (the `sample()` length-1 gotcha). | P2 | S | I-b |
| A1-22 | D1-07 | Unused factor levels in `z` give NaN phi columns. | P3 | S | I-b |
| A1-23 | D1-16 | Each batch reuses the same seed. Results are still order-independent and equal to running each batch alone. | P3 | S | III-a (document) |
| A1-24 | D7-04 | `plotDensity()` smooths each curve with its own bandwidth, so decontaminated peaks look 1.2-1.6x taller. `plot_decontPro.R:61-67` | P3 | S | I-b |
| A1-25 | new | Stale `decontX_<b>_UMAP` entries survive an unbatched re-run; `groupClusters` errors on batched labels. | P3 | S | I-b |
| A1-26 | new | Cosmetic: `runParams` lacks `seed`/`iterLogLik`, typos in messages, `simulateContamination(K = 1)` divides by 0, `retrieveFeatureIndex(NA)` matches an NA row. | P3 | S | I-b |

Not reproduced: D1-10. All supported input classes (dense integer and
double, dgC, dgT, DelayedMatrix over dgC, HDF5) give bit-identical
results.

### Area 2: compiled code

| ID | Finding | Sev | Eff | Phase |
|---|---|---|---|---|
| D2-01 | `_colSumByGroup_numeric` has no NA guard and segfaults; it has no package caller. `matrixSums.c:50-88` | P2 | S | I-d (delete) |
| D2-02 | `decontXInitialize` has no `min(z) < 1` guard: 0 corrupts a header, negatives segfault. `DecontX.cpp:227-238` | P2 | S | I-b |
| D2-03 | The native routine registration is generated, not hand-edited, so the 2026-09 "move it to init.c first" ordering rule is obsolete. | P3 | S | I-d (no action) |
| D2-04 | Dead compiled code: all of `matrixSumsSparse.cpp`, the numeric colSumByGroup, `fastNormPropLog/Sqrt`, `nonzero`. | P2 | S | I-d |
| D2-05/06 | The dead `*ChangeSparse` functions mutate their input in place and don't check dimensions. | P2 | S | I-d (goes with D2-04) |
| D2-07 | `int` offsets overflow at G x C >= 2^31; integer sums overflow silently. `matrixSums.c` | P3 | S | I-d (replaced by `rowsum()`) |
| D2-08 | `decontXEM` divides 0/0 for an empty cluster; R-side guards stop it first. | P3 | S | I-b (with A1-22) |
| D2-09 | `calculateNativeMatrix` exp/log round trip and binary searches: 10x slower than the direct product (max relative diff 7.7e-15). `DecontX.cpp:317-339` | P3 | S | I-e |
| D2-10 | `Eigen::MappedSparseMatrix` is deprecated and removed at Eigen HEAD. The migration compiles clean and is bit-identical. | P2 | S | I-d |
| D2-11 | RcppParallel import lacks a "do not remove" comment. | P3 | S | I-d |
| D2-13 | GPL-3 header in the generated Stan exports, no LICENSE note. | P3 | S | I-d |
| D2-14/15 | `delta` message says "positive integers" but accepts doubles; harmless loop bound. | P3 | S | I-d |
| A2-01 | `load_all()` compiles at -O0, so the EM is 8.7x slower than an installed build. Benchmarks must use an installed build (`dev/bench/` does). | P3 | S | done |
| A2-02 | `//'` roxygen in `matrixNorm.cpp` generates help pages for unexported helpers. | P3 | S | I-d |

### Area 3: tests and coverage

| ID | Finding | Sev | Eff | Phase |
|---|---|---|---|---|
| A3-01 | No stored reference outputs: tests only compare within a run, so nothing can show that a later fix or refactor changed results. | P1 | M | I-a |
| A3-02 | Coverage targets unmet: `R/` 75.1%, 26 of 35 `stop()` paths untested, `plot_decontPro.R` 0%. | P1 | M | I-c |
| A3-03 | Local and codecov coverage denominators differ (local counts the 718-line generated Stan header). | P2 | S | I-c |
| A3-04 | Hollow tests: `expect_s3_class` on unbuilt ggplots, background checks of dimensions only, a decontPro sum identity that a no-op passes, an unasserted legacyInit test. | P2 | M | I-c (background ones in I-a) |
| A3-05 | The two batched fixtures are identical (same default seed), so batch misalignment can't be seen. `test-decon.R:85-86` | P2 | S | I-a |
| A3-06 | Bare `expect_error()` with no pattern (8 sites). | P2 | S | I-a |
| A3-07 | `expect_equivalent` deprecation warning (#55). `test-decon.R:112` | P2 | S | I-c |
| A3-08 | `CreateSeuratObject(dense)` coercion warning in tests. | P3 | S | I-c |
| A3-09 | Five real Stan fits on the same fixture (5.9 of 11.2 s); no Stan-free test of the output processing. | P3 | S | I-c |
| A3-10 | decontPro `ambient_counts` happy path uncovered. | P2 | S | I-a |
| A3-11 | Old scrapper `runUmap` branch can't be covered; scrapper floor should be >= 1.5. `decon.R:995` | P3 | S | I-d |
| A3-12 | Unreachable code in `plot_decontx.R:424-432`. | P3 | S | I-d |
| A3-13 | The `decontXcounts` closures show 0% coverage although called. | P3 | S | I-d |
| A3-14 | The scater-missing `stop()` can't be tested. | P3 | S | I-c |
| A3-15 | Noisy test output (verbose everywhere, Stan chatter). | P3 | S | I-c |
| A3-16 | No snapshots or withr; top-level fixtures; one fixture mutated by a test. | P3 | S | I-a/I-c |
| A3-17 | `test-matrixSums.R` tests dead code; `simulateContamination()` must be pinned before the `rowsum()` swap. | P3 | S | I-a, then I-d |
| A3-18 | C++ guards untested; direct calls cover them in 0.08 s. | P2 | S | I-c |
| A3-19 | No dgC, DelayedArray or HDF5 input test (the A1-06 site). | P2 | S | I-a/I-b |
| A3-20 | `maxIter = 2` in six calls although the fixture converges in 6-22 iterations (about 0.1 s). | P3 | S | I-c |

A draft of about 70 expectations (+2.7 s) took `R/` to 96.6%, so the
targets are reachable. The I-a pinning plan (which behaviour each fix
pins first, and which tests fail until the fix lands) is in the drafts.

### Area 4: code quality

| ID | Finding | Sev | Eff | Phase |
|---|---|---|---|---|
| A4-01 | The decontX argument surface is repeated at about 11 sites with real default drift (`maxIter` 500 vs 200, `varGenes` 5000 vs NULL vs 2000). Not reachable today because internal calls pass every argument. (D1-20) | P2 | M | I-d |
| A4-02 | Validation is scattered and repeated per batch (the root of A1-01). | P2 | M | I-d |
| A4-03 | `.decontX` is 307 lines, `.decontXoneBatch` 190 (D1-21). | P2 | M | I-d |
| A4-04 | `runParams` is used as data, stores the processed batch and lacks `seed`/`iterLogLik`. | P2 | S | I-d |
| A4-05 | Core logic sits in the SCE method (UMAP assembly, background dedup). | P2 | S | I-d/II |
| A4-06 | Three copies of the seed block. | P3 | S | I-d |
| A4-07 | Dead or redundant code: duplicated `decontXEM` calls, inert `stopIter`, `deconMethod`, a duplicate `z=` key. | P3 | S | I-d |
| A4-08 | Repeated `.logMessages` arguments; `cat()` in `celda_functions.R:14`. | P3 | S | I-d |
| A4-09 | Replacing plyr with `match()` is not identical for character `z` with 10 or more labels; the exact replacement is in the drafts. (D3-05) | P2 | S | I-d |
| A4-10 | reshape2 at 5 sites; a faithful base-R helper is identical. | P3 | S | I-d |
| A4-11 | Duplicated plotting code across marker plots and `plotDensity`/`plotBoxByCluster`. | P3 | S | I-d |
| A4-12 | decontPro's three methods are copies; the generic dispatches on `cell_type` too. | P3 | S | II |
| A4-13 | Inconsistent internal names (`.processPlotDecontXMarkerInupt`, snake_case Stan helpers). | P3 | S | I-d |
| A4-14 | Stan helpers lack `@noRd`; decontPro priors documented as "variance" but used as SDs. | P3 | S | I-d / III-a |
| A4-15 | `seq(K)` instead of `seq_len()`, `sum(is.na()) > 0`, and similar lint-level items. | P3 | S | I-d |
| A4-16 | Undeclared `stats`, `grid`, `utils` imports; single-use imports. | P3 | S | I-d |
| A4-17 | After the refactor, move `decontXcounts` and `simulateContamination` into their own files. | P3 | S | I-d |
| A4-18 | Plots read containers directly, so they can't accept other classes. | P2 | M | II |
| A4-19 | Marker plots subset the whole SCE, copying every assay. | P3 | S | II |
| A4-20 | `.process_stan_vb_out` computes the same denominator three times. | P3 | S | I-e |

### Area 5: speed and memory

Measured with `dev/bench/bench-decontx.R` on an installed -O2 build;
simulated data at 10k, 50k and 100k cells, plus PBMC4k and decontPro
on 40 ADTs. Full tables are in the drafts.

| Run | Wall time | Peak RSS |
|---|---|---|
| 10k cells, `z` supplied | 7.6 s | 1.7 GB |
| 50k cells, `z` supplied | 26 s | 2.9 GB |
| 50k cells, 4 batches | 36 s | 3.0 GB |
| 100k cells, `z` supplied | 49 s | 3.7 GB |
| 100k cells, 4 batches | 59 s | 4.6 GB |
| PBMC4k, dgCMatrix | 5.8 s | 1.6 GB |
| PBMC4k, HDF5 | 18.7 s | 2.4 GB |
| decontPro 40 ADT x 5000 cells | 260 s | 1.6 GB |

Peak RSS includes about 0.9 GB for loading the package. Simulated data
converge in 11-16 EM iterations against 90 on PBMC, so EM work is
measured on real data.

| ID | Finding | Sev | Eff | Phase |
|---|---|---|---|---|
| A5-01 | Initialization (PCA, kNN, UMAP) is about 90% of runtime at 50-100k cells and runs even with user `z`. An opt-out when `z` is given is bit-identical for the decontamination. | P2 | S | I-e |
| A5-02 | Opt-in `numThreads`: UMAP 36.7 -> 24.1 s at 100k, bit-identical; PCA 8.4 -> 3.0 s with about 1e-11 differences. | P3 | S | I-e |
| A5-03 | Bounds-checked indexing is about 20% of the EM loop (about 4 s at 100k real-like); unchecked access is bit-identical. | P2 | S | I-e |
| A5-04 | `calculateNativeMatrix` direct form: about 1.7 s and 0.65 GB less peak at 100k, within 1e-14 (D2-09). | P2 | S | I-e |
| A5-05 | Batched assembly makes about 7.5 full copies; one `cbind` and one reorder are identical and cut it to about 2. (D1-27) | P2 | S | I-d |
| A5-06 | HDF5 input: writing the result back to HDF5 is 63% of runtime, and the saved object points at a temporary file. | P2 | S | I-d |
| A5-07 | `anyNA()` on the full counts allocates 217 MB at 100k; checking `@x` is identical. | P3 | S | I-b |
| A5-08 | rstan's CSV in `tempdir()` (about 90 bytes per ADT x cell) is kept until R exits. | P3 | S | I-e |
| A5-09 | SeuratObject costs about 110 MB at load; moving it to Suggests. | P3 | M | II (low priority) |

decontPro: 98% of its time is Stan's ADVI gradient, O(N x M) per
iteration, and the iteration count is unpredictable (2,500 to 12,700;
cap 50,000). #46 is fixed at this scale. No bit-identical speedups.

### Area 6: object coupling and multi-container support

Checked against SeuratObject 5.4.0, anndataR 1.2.2, and ExperimentSubset
1.22.0, then ExperimentSubset was re-checked on 1.23.16 (the PR #54
branch, which includes #53).

| ID | Finding | Sev | Eff | Phase |
|---|---|---|---|---|
| A6-01 | On an ExperimentSubset object, `assayName = "<subset>"` runs on the subset, then fails writing colData; when the sizes divide evenly, contamination is first recycled silently into colData. Validate `assayName` against `assayNames(x)`. | P3 | S | I-b |
| A6-02 | The three ExperimentSubset classes share no class or predicate; `is(x, "SubsetSingleCellExperiment")` misses the spatial and tree classes. Fixed upstream in ExperimentSubset #5; decontX adds no detection code. | P2 | S | II (upstream) |
| A6-03 | `setOldClass("AbstractAnnData")` alone doesn't dispatch S4 methods for anndataR objects. Resolved by decision: an `inherits(x, "AbstractAnnData")` redirect in the ANY methods. | P2 | S | II |
| A6-04 | The contamination plot depends on UMAP column names, which an h5ad round trip drops; `uns` drops names and dimnames. Plots index by position. | P2 | S | I-b (plot), II |
| A6-05 | ExperimentSubset assay names weren't scoped per subset. Fixed upstream by ExperimentSubset #53 (open, after Bioc 3.24). | P2 | S | II (upstream) |
| A6-06 | ExperimentSubset colData row names. Fixed in ExperimentSubset devel (#48). | P3 | S | done upstream |
| A6-07 | `es[i, ]` leaves stale subset indices (ExperimentSubset #4, open). | P3 | S | II (upstream) |
| A6-08 | Seurat details: v3 assays have no `JoinLayers`; reduction keys must be sanitized; `Misc<-` warns on re-run. | P3 | S | II |
| A6-09 | Confirms A1-04 and that decontPro ignores `altExpName`. | P2 | S | I-b/II |
| A6-10 | AnnData views refuse writes; a backed object in mode "r" fails cryptically. | P3 | S | II |
| A6-11 | Per-batch UMAP names are missing from the plan's table. | P3 | S | II (ADR-0005) |
| A6-12 | ExperimentSubset by the direct approach (decision 2026-10-27): read `assay(x, subsetName)` and the subset's colData, call the core, and write back with `setSubsetAssay()`, `colData(x, subsetName =) <-`, `reducedDim(x, , subsetName =) <-` and `metadata(x, subsetName =) <-`. A prototype matched the matrix method for SCE, spatial and tree subsets, siblings, a nested subset, a re-run, and a batched run with `batch` from a main-object column. | P2 | M | II |
| A6-13 | Writing back only the subset's own colData avoids ExperimentSubset's silent drop of overwritten inherited columns (#50). | P2 | S | II |
| A6-14 | If the main object already has `decontXcounts`, `setSubsetAssay()` needs `unlink = TRUE`; it unlinks only that assay. | P2 | S | II |
| A6-15 | The subset's primary assay is the default input, not `"counts"`; reject `assayName = "decontXcounts"`. | P3 | S | II |
| A6-16 | `colData(..., parentColData = TRUE)` returns duplicate names after a subset write (ExperimentSubset #57); read `z`/`batch` columns with that in mind. | P3 | S | II |
| A6-17 | A subset can't see the main object's altExps (ExperimentSubset #58), so decontPro on a subset uses the subset's main assay only (decision). | P3 | S | II |
| A6-18 | The A6-01 silent-recycling mode; covered by the A6-01 validation. | P3 | S | I-b |
| A6-19 | `decontXcounts(es, subsetName =)` forwards to `assay(x, s, subsetAssayName = "decontXcounts")`. | P3 | S | II |

Read and write locations, as verified:

| Class | Input | Decontaminated counts | Per cell | UMAP | Run info |
|---|---|---|---|---|---|
| SCE, SpatialExperiment | `assay(x, assayName)`; decontPro `assay(altExp(x, altExpName), assayName)` | assay `decontXcounts`; decontPro assay `decontPro` in the altExp or main | colData | `decontX_UMAP`, or `decontX_<b>_UMAP` padded with NA | `metadata(x)$decontX`; decontPro in the altExp or main metadata |
| ExperimentSubset | `assay(x, subsetName)` | `setSubsetAssay(x, s, m, "decontXcounts")` | `colData(x, subsetName = s) <-` | `reducedDim(x, , subsetName = s) <-` | `metadata(x, subsetName = s)$decontX` |
| Seurat | Assay5: `JoinLayers` on a copy, then `LayerData` counts; v3: `LayerData` | new assay `decontXcounts` | `meta.data` by `Cells()` | reduction `decontX_UMAP`, key `decontXUMAP_` | `Misc(seu, "decontX")` |
| AnnData (concrete, not views) | `t(X)` or `t(layers[[l]])` | `layers[["decontXcounts"]]` | `obs` | `obsm[["X_decontX_umap"]]` | `uns$decontX` (run parameters and per-batch estimates only) |

Dependencies: anndataR adds nothing beyond decontX's current tree and
goes in Suggests; h5ad needs rhdf5. ExperimentSubset (Suggests) adds 40
packages to decontX 1.11.2's tree. zellkonverter was rejected: it needs
a conda Python environment at runtime and returns an SCE, not an
AnnData object.

### Area 7: docs

| ID | Finding | Sev | Eff | Phase |
|---|---|---|---|---|
| A7-01 | decontX vignette builds `cellTypeMappings` in a hidden chunk, so the shown code fails when pasted; the text says decontX works out the mapping. `decontX.Rmd:163-177` | P1 | S | III-a |
| A7-02 | decontX vignette Seurat section fails three ways: defunct `GetAssayData(slot =)`, `SingleCellExperiment()` not attached, `seuratObject` vs `seuratObj`. | P1 | S | III-a (stopgap PR early in Phase I) |
| A7-03 | decontPro vignette prints about 45 lines of Stan output ("EXPERIMENTAL ALGORITHM", "MAY BE DIVERGING") and explains one of four outputs. | P1 | S | III-a |
| A7-04 | decontPro vignette needs a rewrite (outline in the drafts). | P1 | M | III-a |
| A7-05 | decontPro vignette clusters with Seurat and a dplyr pipe, prints warnings, and downloads 23 MB to use 0.3 MB; cord_blood ships cell-type labels. | P2 | S | III-a |
| A7-06 | decontX vignette runs the full EM twice. | P2 | S | III-a |
| A7-07 | No FAQ in either vignette. | P2 | M | III-a |
| A7-08 | Contamination called a "percentage" but is a proportion (0-1). `decontX.Rmd:118`, `decon.R:105` | P2 | S | III-a |
| A7-09 | Nested `lapply` closures in vignette code. | P2 | S | III-a |
| A7-10 | `simulateContamination()` docs: wrong `delta` default, misspelled return names. | P2 | S | III-a |
| A7-11 | decontPro roxygen: priors called "variance", `ambient_counts` described as N x M (only row sums are used), `cell_type` 1..K rule missing. | P2 | S | III-a |
| A7-12 | Neither vignette nor the package help page cites the papers; the help page cites only RStan. | P2 | S | III-a |
| A7-13 | "Used only when z is not provided" is false. `decon.R:67-90` | P2 | S | I-b |
| A7-14 | `AGENTS.md` gives the DecontPro year as 2023; it is 2024 (NAR 52(1):e4). | P3 | S | III-a |
| A7-15 | `decontXcounts` has no example. | P3 | S | III-a |
| A7-16 | Internal helpers have man pages. | P3 | S | I-d |
| A7-17 | `plotDensity`/`plotBoxByCluster` append ".pdf" to `file` and document the wrong return. | P3 | S | III-a / I-b |
| A7-18 | `retrieveFeatureIndex` roxygen gaps. | P3 | S | III-a |
| A7-19 | Unlabelled chunks, long lines, typos, the decontPro `IndexEntry`. | P3 | S | III-a |
| A7-20 | README lacks a citation section and usage snippet. | P3 | S | III-a |
| A7-21 | DESCRIPTION Title not title case; Description lacks `<doi:>` links. | P3 | S | III-a |
| A7-22 | `z` is duplicated inside `estimates`. | P3 | S | I-d |
| A7-23 | decontPro example uses a 100 x 10 random matrix with no seed; 10 x 8 runs in under a second. | P3 | S | III-a |

Build times (local, cached data): decontX vignette 27 s (two EM runs),
decontPro vignette 27 s (12.6 s in one fit). cord_blood has no empty
droplets.

### Area 8: targeted methods check (issues)

| ID | Finding | Sev | Eff | Phase |
|---|---|---|---|---|
| A8-01 | decontX refuses a single cluster even with `background`, where it is a valid model (correlation with truth 0.9985); the error doesn't name the batch. Allow one cluster per batch that has background. | P2 | S | I-b |
| A8-02 | The k-means fallback (when dbscan finds one cluster) runs on a dense all-genes matrix, the likely cause of #43's 1.4 GiB warning. Run it on the PCA components; changes results on that path only (NEWS); `legacyInit` untouched. | P2 | S | I-b |
| A8-03 | The 2026-09 doublet advice ("run decontX first") is wrong. Simulation: a doublet clustered with one parent loses the other parent's markers, and scDblFinder's AUC on decontaminated counts drops from 1.00 to 0.60. | P2 | S | III-a |
| A8-04 | The drafted #43 workaround `z = rep(1, n)` errors until A8-01 lands. | P2 | S | I-b |
| A8-05 | The #46 fix isn't on RELEASE_3_23; it ships in 3.24. No backport. | P3 | - | decided |
| A8-06 | The 2026-09 D7-03 line references are stale. | P3 | S | done here |

Issue triage: #40 reply and close (FAQ, plus the A1-24 bandwidth fix);
#43 reply (A1-02, A8-01); #44 FAQ (rewrite per the doublet decision);
#45 FAQ; #46 close (fixed in 1.11.1); #39 close (a Bioc 3.21
SummarizedExperiment bug); #56 close (fixed by #63); #55 and #57 stay
open as trackers. #42 (Xenium) is out of scope: the maintainer plans a
new spatial model after this program.

## Maintainer decisions

All decisions, with dates, are in the drafts under each area. The ones
that change behaviour or the plan:

**Results and API (each gets a NEWS entry)**
- Dense input keeps returning a dgCMatrix; the dead coercion is deleted
  and the behaviour documented (A1-14).
- User `z` comes back as a factor with the original labels. Internal
  codes keep today's order so results stay bit-identical; phi/eta
  columns are named by label; `as.integer(result$z)` may change (A1-09,
  A4-09). This lands together with the A1-07 plot fix.
- Initialization still runs by default when `z` is given (unchanged
  output), plus a new opt-out argument; `plotDecontXContamination()`
  accepts any existing reducedDim (A1-16, A5-01).
- HDF5 input returns an in-memory sparse result wrapped as a
  DelayedMatrix; users save to disk explicitly. HDF5Array isn't added
  to Suggests (A5-06).
- Opt-in `numThreads` (default 1); PCA differs by about 1e-11 when it is
  above 1 (A5-02).
- Allow a single cluster with `background` (A8-01); the k-means
  fallback runs on PCA components (A8-02).
- `calculateNativeMatrix` rewrite allowed within about 1e-14 (D2-09).
- Per-batch seed reuse is documented, not changed (A1-23).
- `legacyInit`: no change; scater's deprecation warnings show. An ADR
  follows only when scater makes those functions defunct.
- decontPro: no model or VB-default changes. `iter`, `seed`,
  `output_samples` and a quiet option come with the decontPro-arguments
  ADR; unknown arguments in `...` become an error.

**Code structure**
- One params constructor with no defaults, full method signatures, and a
  test that the methods' formals agree (pattern A, A4-01).
- plyr replaced exactly (internal order kept), reshape2 by a faithful
  helper after the I-b plot fixes; a pure file-move PR after I-d.
- Dead compiled code deleted; the remaining `.colSumByGroup` becomes
  base `rowsum()`, dropping `matrixSums.c` and `R/matrixSums.R`
  (coordinate with celda, which has copies). Generated registration
  stays; `.Call()` uses symbols.
- `cat()` fixed in decontX, and an issue raised in celda.

**Tests**
- `.covrignore` for generated Stan code; coverage targets measured on
  the codecov denominator.
- Committed RDS reference outputs (tolerance 1e-10, user-`z` paths) plus
  a `dev/` script comparing two installed commits with `identical()`
  (no-`z`, legacyInit and decontPro, run locally).
- scrapper >= 1.5 and delete the old-API branch; a `.hasScater()` helper;
  Stan silenced in tests only; regexp or class checks first, snapshots
  for wording; the legacyInit test renamed.

**Multi-container (ADR-0005)**
- ExperimentSubset: the direct approach (A6-12); `decontXcounts(es,
  subsetName =)`; run info in the subset's metadata; decontPro uses the
  subset's main assay only. No ExperimentSubset detection code in
  decontX until ExperimentSubset #5 lands; no workarounds for upstream
  bugs. ExperimentSubset support comes last in Phase II, with a Suggests
  floor on the release that has #53.
- AnnData: a thin anndataR wrapper in Suggests, same object in and out,
  modified in place. It extracts the matrix, resolves `z`/`batch`, calls
  the core (which clusters when `z` is NULL), and writes back to
  `layers`/`obs`/`obsm`/`uns`. Dispatch by an `inherits()` redirect in
  the ANY methods, not `setOldClass`. anndataR's own error for views is
  allowed through. The h5ad -> SCE -> h5ad route is documented as an
  alternative.
- Seurat: a separate `decontXcounts` assay; per-batch reductions hold
  only that batch's cells.
- decontPro on SCE writes its output into the altExp it read.
- decontPro arguments: camelCase is primary; the snake_case names stay
  as permanent aliases (error if both are given). The empty-droplet
  argument becomes `background` (matching decontX), with `ambient_counts`
  a permanent alias. Not `backgroundCounts`/`ambientCounts`, which clash
  with the output elements of those names.

**Docs**
- decontX vignette keeps the full PBMC4k; the second EM run switches to
  the init opt-out once it exists (`eval = FALSE` until then).
- decontPro vignette uses cord_blood's bundled cell-type labels, drops
  Seurat clustering and dplyr (dplyr leaves Suggests), subsamples 500
  cells, and hides Stan output with `results = 'hide'` until the quiet
  option exists.
- No empty-droplet demo in the decontPro vignette or article: decontPro
  doesn't run on empty droplets. `background` gets a short example
  like decontX's (`eval = FALSE`).
- The decontX Seurat section is fixed now, in a small PR early in
  Phase I; III-b replaces it with the native method.
- DecontPro is cited as 2024 everywhere, including `AGENTS.md`.
- Doublets (#44): run decontX and the doublet detectors on the same
  original counts. It is fine when both flag the same cells, and high
  contamination can also be evidence of a doublet.
- Drafted issue replies (#40, #43, #46, #39, #56) are posted once this
  audit is approved.

**Environment**
- Phase 0 ran on the local Bioc 3.23 library. Before Phase I, a separate
  devel library through an untracked repo-root `.Renviron`
  (`R_LIBS_USER`), not renv.

## Open questions (answered 2026-10-27, except 4)

1. **Batch prefix for user `z` in batched runs.** Answered 2026-10-27:
   clusters decontX estimates (no `z`) keep the `"<batch>-<k>"` prefix;
   a user-supplied `z` comes back as the user's labels, unprefixed
   (A1-09). NEWS.
2. **decontPro generic signature.** Answered 2026-10-27: yes. Dispatch
   on the first argument only (`setGeneric(signature =
   "filtered_counts")`), so a call that uses the `cell_type` alias
   doesn't dispatch on a missing `cellType` (ADR-0006); and declare
   `stats`, `grid` and `utils` in Imports (A4-16). Neither changes
   behaviour.
3. **R dependency floor** (#57). Answered 2026-10-27: raising
   `R (>= 4.3.0)` is allowed when needed, for example to clear the
   BiocCheck NOTE or for a dependency; it goes in the dependencies PR.
4. **#44 FAQ and reply text.** The draft in the review is superseded by
   the doublet decision; the rewritten text goes to the maintainer before
   it is posted.

## Changes to the approved plan

Recorded in the plan as amendments (2026-10-27):
- ExperimentSubset: the direct approach replaces the plan's
  `is(x, "SubsetSingleCellExperiment")` branch; no detection code until
  ExperimentSubset #5; the dependency is on ExperimentSubset #53 only.
- AnnData: dispatch by an `inherits()` redirect, not
  `setOldClass("AbstractAnnData")` (which doesn't dispatch); the method
  calls the core directly.
- decontPro: the `ambient_counts` empty-droplet workflow leaves the
  vignette rewrite and the full-data article; `background` gets a short
  example; the argument is renamed `background`.
- HDF5 input returns an in-memory DelayedMatrix (A5-06), a result-class
  change the plan didn't foresee.

## Proposed work

Nothing below is filed until this audit is approved. Each line is one
issue and, normally, one PR, in the order of `dev/ROADMAP.md`.

**Issues to open**
1. Fix the decontX vignette Seurat section (A7-02).
2. Regression net: reference outputs and the cross-commit gate (A3-01,
   A3-05, A3-06, A3-10, A3-16, A3-17, A3-19).
3. Validate `z`, `batch` and scalar parameters (A1-01, A1-10, A1-17,
   A1-18, A1-19, A6-01, A7-13, A5-07).
4. Match `background` rows and check `bgBatch` (A1-02, A1-03, A1-12; #43).
5. Allow one cluster with background; k-means fallback on PCA (A8-01,
   A8-02; #43).
6. Keep user cluster labels, and fix the marker plots (A1-09, A1-07,
   A1-20, A1-22, D2-08).
7. Small crashes: one-cell batches, DelayedMatrix seed, input classes,
   C guard (A1-06, A1-08, A1-14, A1-15, D2-02).
8. decontPro input checks: Seurat split layers, `cell_type`, `...`
   (A1-04, A1-11, A1-13).
9. Plot fixes: shared `plotDensity` bandwidth, stale per-batch UMAPs,
   `file` argument (A1-24, A1-25, A7-17).
10. Cosmetic fixes (A1-21, A1-26).
11. Coverage to target (A3-02 to A3-20).
12. Container-agnostic core (ADR-0004; A4-01 to A4-08, A4-17, A5-05,
    A5-06, A7-22).
13. Drop plyr and reshape2 (A4-09, A4-10).
14. Delete dead compiled code; `rowsum()`; `MappedSparseMatrix`
    migration (D2-01, D2-04 to D2-07, D2-10, D2-11, D2-13 to D2-15,
    A2-02, A7-16).
15. Init opt-out when `z` is given; `numThreads` (A5-01, A5-02).
16. EM and native-matrix speedups (A5-03, A5-04).
17. decontPro arguments: camelCase, aliases, `background`, `iter`,
    `seed`, quiet, temp-file cleanup (decontPro-arguments ADR; A5-08,
    A4-12, A4-14).
18. Seurat methods; AnnData methods; ExperimentSubset methods (ADR-0005;
    A1-05, A4-18, A4-19, A6-*).
19. decontPro vignette rewrite and decontX vignette fixes with FAQ (A7-*,
    A8-03).
20. pkgdown site (pkgdown ADR).

**Issue replies to post** (drafts in the review, area 8): #40 (reply,
close), #43 (reply), #46 (close), #39 (close), #56 (close). #44 waits
for the rewritten FAQ.

**ADR stubs**
- ADR-0004: container-agnostic core (core contract and return shape,
  params constructor, decomposition, file moves).
- ADR-0005: multi-container methods (Seurat, AnnData, ExperimentSubset,
  decontPro returning the object, write locations, per-batch UMAP names,
  dependencies).
- ADR-0006: decontPro arguments (camelCase with aliases, `background`,
  `iter`/`seed`/`output_samples`, quiet option, error on unknown
  arguments, generic signature).
- ADR-0007: pkgdown site (layout, release-only rebuilds, deploy trigger).
