# Comprehensive engineering review — decontX

- **Date:** 2026-09-21
- **Version reviewed:** 1.11.1 (devel branch, commit 92fb6b7)
- **Scope:** engineering health (code, compiled code, dependencies, tests,
  CI, docs, metadata) plus targeted method-level checks on issues #46 and
  #40. No full statistical-methods audit (per scope decision).
- **Method:** local baseline (`make lint/test/check/bioccheck` + covr),
  direct CI-log inspection via `gh`, and seven review passes (one per
  dimension) with runtime verification of claims — every P0/P1 finding
  was either executed live or re-verified against source before landing
  here. This audit changes no package code.
- **Companion document:** `dev/ROADMAP.md` (rewritten from this audit's
  prioritization; the roadmap is the actionable view, this is the
  evidence).

## Executive summary

The package's numerical core is sound and its compiled hot path
(`src/DecontX.cpp`) is well-guarded, but the review found:

1. **Devel CI is fully red.** All four R-CMD-check platforms fail:
   three because decontX's examples trigger celda's new
   `'normalizeCounts' is deprecated` warning, one (macOS) from a runner
   toolchain issue (Rtsne/OpenMP). The lint job also fails (246-lint
   backlog vs. blocking `LINTR_ERROR_ON_LINT`). Local `make check` is
   clean only because the local celda (1.24.0) predates the deprecation.
2. **One confirmed P0 correctness bug:** numeric `batch` labels are used
   as *positional* list indices, silently attributing per-batch estimates
   to the wrong batch (or `NULL`). Verified live with `batch = c(2,1)`
   (silent swap) and `c(3,4)` (NULL + NA-padded names).
3. **Two P0 dependency breaks:** `decontPro()` on Seurat objects errors
   unconditionally (`GetAssayData(slot=)` is defunct in SeuratObject
   5.x; untested, so CI never saw it), and the celda dependency — 40+
   transitive packages — exists for one live call to a function decontX
   already ships an identical C++ copy of (`fastNormProp`).
4. **Issue #46 root cause confirmed with matching math:** rstan::vb's
   default 1000 output draws (two N×M parameters = 99.75% of scalars)
   produce the reporter's ~215 GB CSV / ~317 GB RAM; decontPro uses only
   the posterior-mean row, and `output_samples = 10` is bit-identical at
   ~90× less memory.
5. **Test coverage (65%) materially overstates the safety net:** the two
   heaviest test blocks contain zero assertions ("empty test" skips),
   the only multi-batch fixture is computed and discarded, and nothing
   compares `decontX()` output to the simulator's ground truth.
6. Substantial but cheap wins in docs/metadata: no CITATION for a
   two-paper methods package, NEWS missing six release series (including
   an undocumented breaking rename of `decontPro()`'s first argument in
   1.4.1), no URL/BugReports, a README.md/README.Rmd divergence landmine,
   and no FAQ although 4 of 5 open usage issues are answerable from
   existing behavior.

Severity: P0 = user-facing wrong results/breakage, P1 = real risk or
major maintainability, P2 = worthwhile, P3 = nice-to-have.
Effort: S < half day, M = 1–3 days, L = larger. "ADR" = structural per
AGENTS.md (dependency changes, API/dispatch redesign) — issue first.

## Baseline (2026-09-21, macOS arm64, R 4.5.1)

| Check | Result |
|---|---|
| `make check` (rcmdcheck, --no-manual) | **OK — 0 errors / 0 warnings / 0 notes** (3m35s) |
| `make test` | 31 PASS / 0 FAIL / **15 WARN / 2 SKIP** (skips = empty tests) |
| `make lint` | **246 lints** (commented_code 82, indentation 63, return 32, …) |
| `make bioccheck` (tarball) | 0 ERR / 1 WARN (devel-version false positive) / **10 NOTES** |
| Coverage (covr) | **65.07%** total |
| CI (devel @ 92fb6b7) | **R-CMD-check 0/4 green; lint failing;** coverage + BiocCheck green |

Coverage per file: plot_decontPro.R **0%**, matrixNorm.cpp **0%**,
matrixSumsSparse.cpp **0%**, matrixSums.c 18%, celda_functions.R 30%,
decontPro.R 62%, matrixSums.R 64%, decon.R 75%, plot_decontx.R 81%,
DecontX.cpp 86%, stan_helpers.R 100%.

CI failure causes (verified in run 35633649184 and 35623944284 logs):
- windows-release, ubuntu-release, ubuntu-devel: examples emit
  `'normalizeCounts' is deprecated` (celda) 4× → check WARNING → fail.
- macos-release: transitive dep `Rtsne.so` fails `dlopen` on the GH
  runner (missing OpenMP symbol `___kmpc_dispatch_init_4u`) —
  runner/toolchain issue, not decontX code.
- lint: the 246-lint backlog vs `LINTR_ERROR_ON_LINT: true`.
- ubuntu-devel additionally reports examples > 5s (decontPro 12.2s,
  decontX 9.5s) — five man pages share one example that runs a full EM.

Also observed during baseline: covr/rstantools **regenerated**
`src/stanExports_shrinkage.cc` with `boost::ecuyer1988` where the
committed file has the older `boost::random::ecuyer1988` — the committed
generated Stan exports are stale relative to current rstantools/rstan
(restored via git; regenerate deliberately as a roadmap item).

## Dimension 1 — R code structure & correctness

All bugs below were runtime-verified during review.

| ID | Sev | Eff | Finding |
|---|---|---|---|
| D1-01 | **P0** | S | Numeric `batch` labels: `resBatch[[bat]]` (decon.R:569) indexes positionally; `names(resBatch) <- batchIndex` (:589) pads with NA. `batch=c(2,1)` → silent swap; `c(3,4)` → NULL estimates. Fix: `as.character(batch)` after :438 + regression test. |
| D1-02 | P1 | S | Missing `drop=FALSE` (decon.R:472, :562): 1-cell batch → "subscript out of bounds" crash. |
| D1-03 | P1 | S | `.checkBackground` (:1486-1529) never checks background gene rows/order vs `x`; reordered background silently yields wrong estimates; wrong nrow surfaces an internal C++ message. Direct cause of issue #43. |
| D1-04 | P1 | S | `maxIter`/`iterLogLik`/`convergence`/`seed` unvalidated → cryptic internal errors (`maxIter=0` → "attempt to set an attribute on NULL"). Fix: one `.checkDecontXParams()`. |
| D1-05 | P1 | S | `.processvarGenes` (:1273) checks value before length → "condition has length > 1"; same in `.processdbscanEps` (:1285). |
| D1-06 | P1 | S | User cluster labels returned as integer codes (logical `returnZ` at :444 + factor assign :584-587); NA in `z` silently becomes its own cluster (:1005-1025). |
| D1-07 | P2 | S | Factor `z` with gap levels → NaN `phi` columns (C++ sizes by `max(z)`). Fix: droplevels + recode 1..K. |
| D1-08..12 | P2 | S | `.checkBackground` bgBatch check only under duplicate barcodes; batch length unvalidated vs nC; square zero-Matrix() returns `ddiMatrix`; `.checkDelta` accepts `c(0,0)`; negative counts accepted; `.checkCountsDecon` message unrelated to condition + `sum(is.na())` materializes DelayedArray. |
| D1-13 | P2 | M | decontPro: `cell_type` never validated (1-based contiguity assumed at stan_helpers.R:74); hardcoded `seed=12345`, `iter=50000`; `as.matrix()` densify unguarded. Fix incl. new args → small ADR (merges with D7-01). |
| D1-14 | P2 | S | decontPro methods return invisibly (end with assignment). |
| D1-15 | P2 | S | PCA+UMAP computed unconditionally even with `z` supplied; docs claim otherwise (see D6-04 — decide code vs docs). |
| D1-16 | P2 | S | `with_seed` reuses the *same* seed per batch → identical RNG streams across batches. |
| D1-20 | P1 | M | **ADR:** argument surface duplicated across ~11 sites with real default drift: `maxIter` 500 (methods) vs 200 (internals), `convergence` 0.001 vs 0.01 (:666), `varGenes` 5000/NULL/2000. Fix: one params constructor + delegating methods. |
| D1-21 | P1 | M | `.decontX` (306 lines) / `.decontXoneBatch` (199) monoliths; extract assembly/coercion/seed/EM-loop helpers (same-file, not structural). |
| D1-22..24 | P2 | S | Dead code: 95 commented-out lines (:1171-1265); 5 unreachable internals (+ dead locals, always-true `if`, inert `stopIter` guard, duplicate `z=` list key); `R/package_skeleton_test.R`. Delete all. |
| D1-25 | P2 | S | Duplicate `useDynLib` (bare tags in matrixSums.R vs `.registration=TRUE`); string `.Call()` instead of registered symbols. |
| D1-26..27 | P2 | S | 3× duplicated seed blocks (one helper); assembly loop copies full sparse matrix up to 3×/batch. |
| D1-29..30 | P3 | S | Misspelled `.processPlotDecontXMarkerInupt`; stray `#' NULL`; decontPro generic lacks `@export`; vendored celda helpers lack source-version header (coordinate upstream, don't diverge). |

## Dimension 2 — Compiled code (src/)

| ID | Sev | Eff | Finding |
|---|---|---|---|
| D2-01 | P1 | S | `_colSumByGroup_numeric` lacks the NA-group guard its integer twin has → **reproducible segfault** (verified). Same gap in `_rowSumByGroup_numeric`. Currently reachable only from dead code, but one factor(levels=seq(K)) away from live. |
| D2-02 | P1 | S | `decontXInitialize` is the only export missing the `min(z) < 1` guard → out-of-range `z` silently corrupts the R object header (verified). |
| D2-03 | P1 | S | matrixSums registration is **hand-edited into generated `RcppExports.cpp`**; `compileAttributes()` rewrites it; git history shows it mangled twice (573ed54, 1c9ba96). Fix: hand-owned `src/init.c` (or convert to Rcpp after deletions). **Must land before any change that regenerates RcppExports** (e.g. D3-15). |
| D2-04 | P2 | S | 6 of 8 matrixSums.c entry points unreachable (not in registration table; "not in load table" verified) ≈ 250/382 lines; 3 `*Sparse` variants also dead. Commented-out tests calling them can never be re-enabled. Delete. |
| D2-05/06 | P2 | S | Latent (dead-code) bugs: `*Change*` fns mutate caller SEXP in place; ChangeSparse validates only one `px` dimension (unchecked heap write) + wrong error message. Deletion retires both. |
| D2-07 | P2 | S | matrixSums.c `int` offset arithmetic overflows on >2^31-element dense matrices (reachable). Fix: `R_xlen_t`. |
| D2-08 | P2 | M | `decontXEM` divides by possibly-zero cluster colsums → silent NaN (empty cluster); `decontXInitialize` immune (pseudocount pre-fill) — inconsistent hardening. |
| D2-09 | P2 | S | `calculateNativeMatrix` re-introduces the exp(log+log) round-trip `decontXEM` explicitly removed + O(log nnz) `coeffRef` in inner loop + compiler-flagged dead var. |
| D2-10/D3-15 | P1 | M | `Eigen::MappedSparseMatrix` at 8+8 hand-written sites (+16 generated mirrors): doxygen-deprecated, zero warnings today, but removal in Eigen 3.5/4 = total build failure on all platforms at once. Mechanical `Map<SparseMatrix<double>>` swap + `compileAttributes()`; **after** D2-03. |
| D2-11 | P3 | — | Lead correction: **RcppParallel is NOT removable** — rstantools/StanHeaders TBB scaffolding in src/Makevars. Add a do-not-remove comment. |
| D2-13 | P3 | S | GPL-3 header in generated `stanExports_shrinkage.h` vs MIT license: standard rstantools template output; do not edit generated file; note in LICENSE/inst/NOTICE, optionally raise upstream. |
| D2-14/15 | P3 | S | `delta` message says "positive integers" but accepts 0/doubles; eta loop iterates over `new_phi.ncol()`; dead `nonzero()` O(n²) push_back — delete with matrixNorm cleanup (keep `fastNormProp`, see D3-01). |

Positives: `decontXEM`/`decontXLogLik`/`calculateNativeMatrix` each run
6–8 precondition checks; PROTECT/UNPROTECT balanced everywhere;
Makevars portable.

## Dimension 3 — Dependencies & retired APIs

| ID | Sev | Eff | Finding |
|---|---|---|---|
| D3-01 | **P0** | S | **ADR:** celda (40+ transitive pkgs) has ONE live call site (decon.R:1459, `normalize="proportion"`); the other 4 are in dead functions. decontX's own `src/matrixNorm.cpp fastNormProp` is an exact replacement (celda calls the same routine internally). Fix: internal `.normalizeProportion()` + delete dead fns + drop celda. **Resolves the CI failure at the root.** |
| D3-02 | P1 | S | celda exports 6 names identical to decontX's exports; hard dep guarantees co-installation → attach-order roulette. Dropping celda removes the guarantee; file upstream celda issue (re-export, don't duplicate). |
| D3-03 | **P0** | S | `decontPro,Seurat` is broken on current installs: `GetAssayData(slot=)` **defunct** in SeuratObject 5.x (verified error). Fix: `SeuratObject::LayerData(layer="counts")` + guarded regression test. |
| D3-04 | P1 | M | **ADR:** Seurat (45 transitive) → SeuratObject (13): everything used lives in SeuratObject. Land with D3-13 (removing Seurat unhides the broken `\link[uwot]`). |
| D3-05/06 | P2 | S/M | **ADR:** drop plyr (2 sites ≡ `match(z, unique(z))`) and reshape2 (5 melt sites → 2 base-R helpers; regression risk in factor-vs-int columns — gate on new plotting tests). |
| D3-07 | P2 | S | MCMCprecision's only live use is a C++ callback (`namespace_env`, DecontX.cpp:61-63) — invisible to static analysis; move `@importFrom` to package doc with do-not-remove comment. (Also: callback crosses C++→R once per EM iteration — perf note.) |
| D3-08/09 | P2 | S | **ADR:** scran Suggest is dead (only in commented-out block); HDF5Array used undeclared (decon.R:612-614) — prefer deleting the redundant special case (generic canCoerce fallback handles it). |
| D3-10/11 | P2 | S | DESCRIPTION: add URL/BugReports; expand biocViews (add at least RNASeq, Preprocessing, Normalization, QualityControl, Proteomics — all vocabulary-verified). |
| D3-12 | P1 | S | `aes_string` 5 sites (deprecation warnings fire in tests today; 3 sites are literal names → plain `aes()`); also `trans=` → `transform=` (2 sites, silent today). |
| D3-13 | P2 | S | `man/decontX.Rd` `\link[uwot]{umap}` stale (undeclared pkg + factually wrong post-scrapper) → `\link[scrapper]{runUmap}`; typo "dimenions". |
| D3-16 | P3 | S | **ADR (minor):** Matrix >=1.5.3 floor unenforceable (drop); verify scrapper >=1.2.0 floor covers `mean.filter`; R >= 4.3.0 → 4.5.0 on devel. |
| D3-17..19 | P3 | S | withr referenced two ways (normalize); `importFrom(rstan,sampling)` is unused rstantools boilerplate (comment, don't delete); .Rbuildignore gaps (BiocCheck dir, tarballs, vignette html). Good news (verified, don't "pre-fix"): ggplot2 4.0 S7 is a non-event here; `as(., "dgCMatrix")` not deprecated; no `:::` anywhere; all 21 Imports have ≥1 live use. |

**Target dependency state:** Imports 21 → 17 (drop celda, plyr,
reshape2; Seurat→SeuratObject), Suggests −scran (+HDF5Array only if the
special case is kept). Roughly 60–70 transitive packages removed.
Proposed ADR grouping: ADR-D1 celda removal; ADR-D2 Seurat→SeuratObject;
ADR-D3 plyr+reshape2; ADR-D4 optional backends (scran/HDF5Array/dead
init block); ADR-D5 DESCRIPTION metadata/floors.

## Dimension 4 — Tests & coverage

| ID | Sev | Eff | Finding |
|---|---|---|---|
| D4-01 | P1 | M | The two heaviest blocks (test-decon.R:98, :119) execute 6 decontX + 7 plot calls with **zero assertions** — reported as "empty test" skips while silently providing most SCE/background/plot coverage. |
| D4-02 | P1 | S | `modelDecontXoneBatch`/`batchDecontX` fixtures computed and never used; the latter is the suite's ONLY multi-batch call and uses the one numeric labelling that hides D1-01. |
| D4-03 | P1 | M | No behavioral oracle: simulator ground truth never compared to output; decontX returning its input unchanged would pass the suite. Fix: recovery test (correlation with truth; decontaminated closer to native than observed; EM log-lik non-decreasing). |
| D4-04 | P2 | S | `maxIter = 2` in 6/7 calls saves 0.08s (default converges in 9 iters / 0.25s — measured) while making convergence/log-lik branches unreachable. |
| D4-05 | P1 | M | decontPro: 1 test, 1 assertion, **14.8s (98% of suite wall clock)**, Pareto k=12.4 (meaningless fit). 10×8 fixture = 1.9s, k=3.3. Untested: ambient path, both guards, SCE + Seurat methods, parameters, visibility. Use `local_mocked_bindings(.call_stan_vb=…)` for Stan-free dispatch tests. |
| D4-06..08 | P2 | S/M | No helper/setup files; test-matrixSums.R mutates a shared fixture between blocks; zero withr in tests; `expect_error()` without regexp (the `.checkDelta` tests pass on ANY error); ~20+ user-facing messages unpinned; `_snaps/` empty. Policy: `expect_snapshot(error=TRUE)` for message wording; never snapshot numeric estimates. |
| D4-09..13 | P2 | S/M | Highest-value untested surfaces: `retrieveFeatureIndex` (exported; examples cover more than tests), `colSumByGroupSparse` (live, 0%), `.calculateDecontXBarplotPercent` (pure arithmetic + 2 latent hazards), plot_decontPro.R (0%, exported, vignette-featured), plot_decontx `groupClusters` machinery. |
| D4-14 | P3 | S | Plot-test policy: **no vdiffr** (Bioc-builder fragility); assert structure/data via `layer_data`/`ggplot_build`; write policy into AGENTS.md. |
| D4-15 | P2 | S | Coverage CI: no codecov.yml (default target:auto fails on stochastic jitter; denominator counts generated files); `needs: coverage` resolves to nonexistent field; no `make coverage`; upload-artifact@main unpinned; checkout@v3; R-CMD-check lacks `_R_CHECK_FORCE_SUGGESTS_=false` (disagrees with make check). Proposed codecov.yml in review transcript (target 80%, threshold 1%, ignore generated). |
| D4-16 | — | — | The commented-out matrixSums tests are permanently unrunnable (symbols unregistered) — delete with D2-04, never restore. |

**A full 18-file test plan** (fixtures, assertions, runtime class, which
finding each file guards) was produced during review — see the
"Test plan" appendix pointer below. Projection: dead-code deletion alone
(~300+ uncovered lines) moves 65% → high-70s; the plan makes the number
honest and **halves** suite wall clock (16s → ~8s) via the decontPro
fixture fix. Effort ~20–24h total.

## Dimension 5 — CI & repo hygiene

| ID | Sev | Eff | Finding |
|---|---|---|---|
| D5-01 | P1 | — | Devel CI fully red (see Baseline). Wave 1 must make devel green. |
| D5-02 | P2 | S | R-CMD-check.yaml: checkout@v3 (Node-20 deprecation firing); oldrel-1 commented out — decide policy and delete or restore deliberately; add `_R_CHECK_FORCE_SUGGESTS_=false`. |
| D5-03 | P2 | S | test-coverage.yaml: checkout@v3, upload-artifact@main, no CODECOV_TOKEN/fail_ci_if_error, phantom `needs: coverage` (fix with D4-15). |
| D5-04 | — | — | Positive: bioccheck.yaml is well-designed (Bioc devel container, tarball+GitClone, weekly cron, documented flip-to-required plan). |
| D5-05 | P3 | S | No issue templates (would triage the support-question influx), no dependabot (how checkout@v3 lingered), no CODEOWNERS (optional). |
| D5-06 | P2 | S | `dev/agent-log.md` referenced by dev/AUDIT.md but missing — create stub or amend AUDIT.md; the audit process is unrunnable as written. |
| D5-07 | P3 | S | .Rbuildignore/.gitignore gaps (= D3-18). |

Lint strategy: commented_code (82) largely evaporates with the dead-code
deletions; then indentation (63) + return (32) are mechanical. One
style-only PR after deletions; keep lint blocking thereafter.

## Dimension 6 — Documentation & community

| ID | Sev | Eff | Finding |
|---|---|---|---|
| D6-01 | P1 | S | No inst/CITATION; package Rd cites only RStan; vignettes never cite either paper. Draft CITATION produced (Yang 2020 Genome Biology 21:57; Yin NAR 52(1):e4 — pick 2023 vs 2024 phrasing). |
| D6-02 | P1 | S | NEWS missing 1.0.0–1.10.0 incl. an **undocumented breaking rename** in 1.4.1 (decontPro `object` → `filtered_counts`, shipped as a RELEASE_3_20 patch). Draft backfill produced (ancestry git-verified). Add a NEWS checkbox to the PR template. |
| D6-04 | P1 | S/M | "Used only when z is not provided" is false for varGenes/dbscanEps/legacyInit (= D1-15): init always runs; z-supplied users still pay HVG/PCA/UMAP and need scater under legacyInit. Decide docs-match-code (cheap) vs code-match-docs (better; changes returned object) — ADR-level decision. |
| D6-09 | P1 | S | README.md hand-edited despite "generated from README.Rmd" header — re-knitting reverts the Bioc install instructions and the published-paper link. Port to Rmd + reknit, or delete Rmd (it has no computed chunks). Use the DOI-form NAR link; add Bioc build badge. |
| D6-12 | P1 | M | No FAQ surface; the vignette never mentions `batch` or doublets, yet 4 of 5 open usage issues (#40 #43 #44 #45) are answerable from existing behavior. FAQ drafts produced — **with the #40 answer corrected by dimension 7 (see D7-03): decontaminated counts can never exceed observed; do not ship the earlier "model expectation" framing.** |
| D6-05..08 | P2/P3 | S | decontXcounts.Rd example; 6 internal man pages pollute the index (`@noRd`); decontX `@return` under-documents the 10-slot estimates list + magic `"all_cells"` key; decontPro hardcoded seed/iter undocumented + `@return` doesn't name the three matrices. |
| D6-13 | P2 | S | decontPro vignette: redundant untuned chunk; 3 inconsistent marker sets; **1-based `cell_type` footgun** (Seurat clusters are 0-based; `as.integer(Idents())` is load-bearing and invisible) — document + validate (with D1-13). |
| D6-14/15 | P3 | S | No CODE_OF_CONDUCT (link Bioc CoC); issue templates (with D5-05); package-Rd stub; unanchored `\link{assay}`/`\link{Matrix}`. |
| D6-11 | P2 | M | Vignettes: 2 ExperimentHub downloads + 3 model fits per build per platform (decontX.Rmd runs the EM twice). `.Rmd.orig` precompute warranted (decontPro certainly; decontX arguably). Interim: `eval=FALSE` the second decontX run. |

## Dimension 7 — Method-level checks (issues #46, #40)

**#46 (decontPro memory blowup) — confirmed, fix in code.**
decontPro consumes only `@sim$est` (the ADVI posterior-mean row,
stan_helpers.R:28); the 1000 default output draws are computed, written
to CSV, read back, copied twice, and discarded. `delta` and `background`
are both N×M (99.75% of 22.5M scalars at the reporter's 215×52,222) →
CSV ≈ 225 GB and peak RAM ≈ 360 GB — matching the reporter's observed
215 GB / 317 GB within 5%. Measured: `output_samples = 10` yields
**bit-identical** `sim$est` at ~90× less memory. `pars`/`include`
filtering was tested and does **not** reduce the CSV (filtering happens
after the full read); there is no generated-quantities block to move.
`output_samples = 1` crashes in loo::psis — validate ≥ 2. Residual:
rstan runs a per-column PSIS pass (~1h at that scale) that decontPro
cannot disable — document. Fix: default `output_samples = 10` in
`.call_stan_vb` and expose `output_samples`/`seed`/`iter` as `decontPro`
arguments (small ADR, merges D1-13/D6-08). Deeper scaling (2·N·M latent
parameters) is methods-level — separate tracking issue, explicitly a
non-goal without a methods review.

**#40 (increased decontaminated signals) — expected behavior + plot
artifact; close with corrected explanation.**
`decontaminated_counts = counts × rate/(sum of three rates)` elementwise
(stan_helpers.R:88-100), all rates non-negative by Stan constraints →
**decontaminated ≤ observed is a mathematical guarantee**, and the three
components sum exactly to the observed counts (verified numerically;
same guarantee for decontX via DecontX.cpp:338-339). The reported
"spikes" are (a) a per-group KDE bandwidth artifact in `plotDensity`
(decontamination collapses the negative population's IQR → the
decontaminated curve is smoothed less → +24% taller peak measured with
no value increasing), and (b) for normalized values, a legitimate CLR
effect (smaller geometric-mean denominator) — the desired outcome.
Genuine plot defect worth fixing: `plotDensity` should pass a shared
bandwidth to both curves. Add the invariant regression test
`all(decontaminated_counts <= counts)` + the sum identity.

Issue dispositions: #46 fix in code; #40 close with corrected
explanation + FAQ + plotDensity bandwidth fix; #39 close (stale library;
`BiocManager::valid()`); #43/#44/#45 answered by the FAQ + the
background-row validation fix (D1-03).

## Cross-dimension reconciliations

- `fastNormProp` (matrixNorm.cpp): dimension 2 proposed deleting the
  file as dead; dimension 3 makes `fastNormProp` the celda replacement.
  **Resolution: keep `fastNormProp` (mark internal), delete only
  `fastNormPropLog`/`fastNormPropSqrt`/`nonzero` + orphan man pages.**
- #40 FAQ: dimension 6's draft answer is refuted by dimension 7's code
  walk — ship the D7 version.
- RcppParallel: flagged as possibly removable in early exploration;
  dimension 2 verified it is load-bearing (Stan/TBB). Keep.
- Ordering constraint: D2-03 (move registration out of generated
  RcppExports.cpp) must land **before** anything that reruns
  `compileAttributes()` (the Eigen migration D3-15, and even
  `make docs`, which already silently mangles it).

## Verification

- Every P0/P1 was runtime-verified during review (executed failing code,
  or CI logs read directly) or source-verified line-by-line by the lead
  reviewer; the ledger records which. No P0/P1 was downgraded.
- Baseline commands are re-runnable: `make lint/test/check/bioccheck`
  plus `covr::package_coverage()`.
- This audit changed no package code; `git status` is clean apart from
  two baseline-run artifacts (`decontX.BiocCheck/`,
  `decontX_1.11.1.tar.gz`) noted for removal.

## Appendix

Reusable drafts produced during review — NEWS.md backfill, inst/CITATION,
FAQ sketches (with the #40 answer corrected per D7-03), proposed
codecov.yml, and the condensed 18-file test plan — are preserved in
`dev/audits/2026-09-review-drafts.md`. They get materialized by the
roadmap waves that own them.
