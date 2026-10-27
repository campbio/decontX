# 2026-10 deep review: working drafts (Phase 0)

Raw per-area notes from the Phase 0 review (plan:
`dev/plans/2026-10-18-deep-review-and-multi-container-design.md`), run
2026-10-18 to 2026-10-20 in sequential subagents. They are the input to
`dev/audits/2026-10-deep-review.md`, which hasn't been written yet.

These notes were kept only in a session scratchpad, which was later
cleared. They were rebuilt on 2026-10-26 by replaying the file writes
recorded in that session's transcript, so the text matches what was
saved then. Scripts and outputs the notes mention ("scripts in this
dir", `.rds` files, `area6/lib`) were not recovered and no longer exist.

## Area 7 and area 8 decisions

These were still open when the review session ended. The maintainer
answered them on 2026-10-26; the answers follow the two lists.

Area 7 (docs), recommendations:
1. decontX vignette data: keep the full PBMC4k (about 20 s once the
   second EM is gone); subsample to 2000 cells only if the builders
   complain.
2. Second EM run: switch to the skip-initialization option once it
   exists (Phase I); if the docs PR lands first, `eval = FALSE` for now.
3. decontPro cell types: use cord_blood's bundled labels; drop the
   Seurat clustering and the dplyr pipe; remove dplyr from Suggests.
4. Empty-droplet demo: simulated droplets in the vignette; a real
   raw-data example in a pkgdown article.
5. decontPro subsample: 500 cells (about 6 s per fit, three fits).
6. Stan console output: `results = 'hide'` now; a quiet option comes
   with the decontPro-arguments ADR.
7. Fix the broken Seurat section now as a small stopgap PR early in
   Phase I; III-b replaces it.
8. Citation year 2024 everywhere, and fix `AGENTS.md`.

Area 8 (targeted methods check):
1. A8-01: allow a single cluster when `background` is supplied
   (recommended: yes).
2. A8-02: run the k-means fallback on PCA components, not a dense
   all-genes matrix (recommended: yes; NEWS entry; `legacyInit`
   untouched).
3. A8-05: backport the #46 fix to RELEASE_3_23 (recommended: no; it
   ships in 3.24).
4. When to post the drafted issue replies (#40, #43, #46, #39, #56):
   with the audit's approval, or sooner (#43's reporter is blocked).
   Also confirm the #44 doublet-ordering framing (A8-03).

### Maintainer decisions (2026-10-26)

- Area 7: items 1-3 and 5-8 AGREED as recommended.
- Area 7 item 4: **no empty-droplet demo** is needed ("Don't need
  empty drops"). This drops the `ambient_counts` empty-droplet workflow
  from the decontPro vignette rewrite (plan section III-a, item 5).
- Area 8 items 1-3 AGREED as recommended (A8-01 yes, A8-02 yes, no
  A8-05 backport).
- Area 8 item 4: post the drafted replies with the audit's approval.
- #44 doublets (A8-03): **run decontX and the doublet detectors on the
  same original counts**, not one before the other. They sometimes flag
  the same cells, and that is fine. High contamination can also be
  evidence that a cell is a doublet. This replaces both the 2026-09
  "decontX first" draft and the area 8 "remove doublets before decontX"
  draft FAQ below; the FAQ and the #44 reply are rewritten to match.

Since then, ExperimentSubset has changed a lot (PRs #38-#51 merged;
#53 subset assay scoping and #54 `getSubsetExperiment()` /
`setSubsetExperiment()` open, held until after Bioc 3.24). The area 6
ExperimentSubset findings were made against ExperimentSubset 1.22.0
and need rechecking against that API.



---

<!-- recovered from scratchpad/baseline-summary.md -->

## Baseline at c6fca8e (2026-10-18), R 4.6.1, local Bioc 3.23 (release; 35 pkgs out of date), BiocCheck 1.48.1
- lint: no lints
- test: FAIL 0 | WARN 2 | SKIP 0 | PASS 84 (~33s). Warnings: CreateSeuratObject "Data is of class matrix. Coercing to dgCMatrix" (test-decontPro.R:118 + 1 more)
- coverage 74.21%: plot_decontPro.R 0; matrixSumsSparse.cpp 0; matrixNorm.cpp 17.24; celda_functions.R 30.26; matrixSums.R 64.29;
  stanExports_shrinkage.h 80.5; decon.R 83.66; plot_decontx.R 85.11; DecontX.cpp 85.90; decontPro.R 95.77; matrixSums.c 97.73; stan_helpers.R 100
- check: Status OK
- bioccheck: 1 ERROR (support site HTTP 504 -> env/transient), 1 WARNING (y even in release -> env: local BiocCheck is release 3.23; false positive on devel),
  NOTES: R dep 4.3.0->4.6.0 (#57); no fnd role; missing chunk labels both vignettes (#57); sapply in generated R/stanmodels.R (false positive, generated);
  cat() in R/celda_functions.R:14 (defect); 14 functions >50 lines, .decontX 307, .simulateContaminatedMatrix 114 (#57, D1-21);
  decontXcounts.Rd no examples (D6-05); 51 lines >80 chars (incl. generated RcppExports.R, vignettes); 1217 lines not 4-space (2-space policy, .lintr; false positive);
  mailing list (env)
- devel CI green (#63 merge, 2026-10-09)
- git tree clean after all runs (no stanExports churn this time)
- Env gap to flag: local library is Bioc 3.23, package targets devel 3.24.

---

<!-- recovered from scratchpad/2026-09-status-notes.md -->

## 2026-09 finding status at c6fca8e (agent report, 2026-10-18)
Counts: 23 FIXED (incl. D5-01, confirmed by green devel CI on #63), 11 PARTIAL, 40 OPEN, 1 N/A (D5-04).
FIXED: D1-01, D1-14, D3-01, D3-02(decontX side), D3-03, D3-04, D3-07, D3-10/11, D3-12, D3-13,
 D4-01, D4-02, D4-03, D4-15, D4-16, D5-01, D5-02, D5-03, D5-06, D6-02, D6-09, B2-01, CI-01
PARTIAL: D1-22..24 (left: dup z= key decon.R:585, inert stopIter 737/780, always-true deconMethod if 747);
 D2-04 (3 dead *Sparse variants in matrixSumsSparse.cpp); D2-05/06 (Sparse Change variants mutate px :115, msg :164);
 D4-04 (maxIter=2 everywhere except oracle); D4-05 (left: ambient_counts happy path, parameters content);
 D4-06..08 (fixture mutation test-matrixSums.R:21; bare expect_error test-matrixSums 13-14,24-25, test-decon 122,284,288; no withr; no snaps);
 D4-09..13 (untested: retrieveFeatureIndex, colSumByGroupSparse, .calculateDecontXBarplotPercent values, plot_decontPro.R);
 D6-01 (package Rd/vignettes don't cite papers); D6-14/15 (no CoC/templates, minimal package Rd);
 D7-01 (output_samples/seed/iter not exposed); D7-05 (decontPro <= counts not asserted)
OPEN R: D1-02 (decon.R:479-480 drop=FALSE), D1-03 (#43, decon.R:1268-1315), D1-04 (779,822), D1-05 (1050,1062),
 D1-06 (453,591-594,894-897), D1-07 (894; DecontX.cpp:238), D1-08..12 (1298-1306,449,410,1070,867-873),
 D1-13 (decontPro.R:159-166; stan_helpers.R:20-21), D1-15/D6-04 (715-722), D1-16 (521-523), D7-04 (plot_decontPro.R:61-67)
OPEN C: D2-01 (matrixSums.c:50-86), D2-02 (DecontX.cpp:228-238), D2-03 (RcppExports.cpp:184-203), D2-07 (matrixSums.c:34,41,71,80),
 D2-08 (DecontX.cpp:114-118,123-126), D2-09 (334-339), D2-10/D3-15 (MappedSparseMatrix), D2-11, D2-13, D2-14/15 (matrixNorm.cpp:35-105)
OPEN quality: D1-20, D1-21, D1-25 (dup useDynLib, string .Call matrixSums.R:14-24), D1-26..27 (seed blocks 502-541,1024-1039,1123-1146; copies 563-574), D1-29..30
OPEN deps: D3-05/06 (plyr decon.R:878,895,1182; reshape2 plot files), D3-08/09 (scran Suggest; HDF5Array decon.R:618-620), D3-16 floors, D3-17..19
OPEN CI: D5-05 (templates/dependabot/CODEOWNERS), D5-07 (=D3-18)
OPEN docs: D6-05..08, D6-11 (decontX.Rmd:89,265), D6-12 FAQ, D6-13 (decontPro.Rmd:99), D7-03 (#40 still open)
Other: untracked src/stanExports_shrinkage-*.o.tmp (build artifact, no ignore rule); examples run full EM in 5 man pages.
Open GH issues: 57, 56 (may be fixed by 52de8a0), 55, 46, 45, 44, 43, 42, 40, 39.

---

<!-- recovered from scratchpad/area1/findings.md -->

## Area 1: R code bugs (agent report, 2026-10-18; A1-01 and A1-06 re-verified by lead)

Env: decontX 1.11.2 @ c6fca8e, R 4.6.1, SeuratObject 5.4.0, scrapper 1.6.3. Scripts in this dir; setup.R = load_all;
setup_pro.R = 10x8 ADT fixture. decontX fixture simulateContamination(C=100,G=60,K=3,seed=1). All live-reproduced.

| ID | 09 ID | Title | Sev | Where | Repro | Fix | Eff | Phase |
|---|---|---|---|---|---|---|---|---|
| A1-01 | new | z shorter than ncol (incl. column name "cluster", z=2, z[1:90]) runs silently; NA pad -> own cluster; cor w/ truth -0.02 vs 0.999 | P1 | decon.R:499-501 (length check only at :880) | s08_colname_args.R, s05b_silent.R | check length(z)/length(batch)==ncol in .decontX before batch subset | S | I-b (names: II) |
| A1-02 | D1-03 #43 | background rows never matched to x: reversed rows -> contamination max diff 0.47, cor 0.988->0.801; bg[-1,] or superset -> exact #43 error ("rows in 'eta'..."); non-overlapping names run silently | P1 | decon.R:1268-1315, 767-770 | s02_background_rows.R | match by rownames, error on missing; unnamed -> require equal nrow | S | I-b |
| A1-03 | D1-08 up | bgBatch mismatch silently drops background (uniform eta, cor 0.988->0.746); missing batch -> uniform eta; short recycled; long -> cryptic | P1 | decon.R:480, 1298-1306 | s05_d108_12.R, s05b_silent.R | require length(bgBatch)==ncol(bg) and every batch has bg cols | S | I-b |
| A1-04 | new (design) | decontPro on Seurat v5 split layers: LayerData warns + returns first layer (10x4); length-8 cell_type -> Stan dim mismatch then "no applicable method for '@'... list"; length-4 -> runs on half the cells | P1 | decontPro.R:93,95 | s07_seurat_split.R | JoinLayers in memory (or stop) | S | I-b (full: II) |
| A1-05 | new | decontPro reads default Seurat assay / SCE main "counts": ran on RNA not ADT (sum 383 vs 3265); altExp ADT ignored; assay "raw" errors | P1 | decontPro.R:70,72,93,95 | s07, s09_pro_colname.R | assay/altExp args (or warn) | M | II |
| A1-06 | new | DelayedMatrix with dense in-memory seed crashes after full EM: "the condition has length > 1" (seed class is c("matrix","array")) | P1 | decon.R:617-618 | s11_input_types.R, s11b_delayed.R | inherits()/is(), or delete special case (D3-09) | S | I-b |
| A1-07 | new | plotDecontXMarkerPercentage wrong for numeric labels not 1..K: 0-based -> error 9 vs 6 rows; 10/11/12 -> NA labels; 2/5/9 mislabelled; factor w/ unused levels -> wrong values (Inf) | P1 | plot_decontx.R:225-229, 572-581 | s12_plots.R, s12b_plot_pct_check.R | match by name, droplevels, align totals | S | I-b |
| A1-08 | D1-02 | 1-cell batch or 1-cell complement crashes ("subscript out of bounds"); 1-col bg batch -> coercion error | P1 | decon.R:479,480, **569** | s01_onecell_batch.R, s01b | drop=FALSE x3 + explicit 1-cell error | S | I-b |
| A1-09 | D1-06 | user z labels returned as codes (char -> 1 2 3); factor codes single-batch but labels when batched (inconsistent); char NA own cluster; factor NA -> "negative length vectors" | P1 | decon.R:453, 589-594, 894-897 | s04_zlabels.R | keep labels, reject NA | S | I-b (decision: celda/sctk effect) |
| A1-10 | D1-04/05 | scalar params unvalidated: maxIter 0/-1 cryptic; convergence="x" stops after 1 iter silently; maxIter="a" runs; iterLogLik 0/NA; seed length 2; varGenes/dbscanEps length 2; estimateDelta="yes" -> FALSE; maxIter<iterLogLik -> logLikelihood NULL | P1 | decon.R:779-822, 1050, 1062 | s03_params.R, s17_misc.R | one .checkDecontXParams() | S | I-b (consolidate I-d) |
| A1-11 | D1-13+ | decontPro cell_type/counts unvalidated; Stan failure (vb returns list) -> "no applicable method for '@'"; gap (1,3) "no more scalars"; NA -> "not a matrix"; seed/iter hardcoded | P2 | decontPro.R:159-166; stan_helpers.R:15-30,41 | s06_decontpro_celltype.R | validate/recode cell_type to 1..K; stop if not stanfit | S | I-b (seed/iter ADR w/ D7-01) |
| A1-12 | new | decontPro ambient_counts rows not matched (decontPro #43): reversed accepted, p["ADT1"] 0.671->0.034; real fit max diff 0.93 | P1 | decontPro.R:137-147 | s14_pro_ambient.R, s15_pro_zero_row.R | shared row-matching helper with A1-02 | S | I-b |
| A1-13 | new | decontPro `...` swallows unknown args (seed=1 vs 2 identical; delta_SD typo ignored; batch ignored) | P2 | decontPro.R:69,92,114 | s06/s09 | error on non-empty ... | S | I-b |
| A1-14 | new | dense input returns dgCMatrix (canCoerce(..., c("matrix","array")) errors in silent try) | P2 | decon.R:630-638 | s01, s11 | class(counts)[1]; dense-out is a decision | S | I-b |
| A1-15 | new | as(x,"dgCMatrix"): dgRMatrix error; dgTMatrix Matrix deprecation warning; data.frame/SE/matrix x w/ SCE bg cryptic | P2 | decon.R:491,495 | s11 | as(., "CsparseMatrix") + input class check | S | I-b |
| A1-16 | D1-15/D6-04 | init (HVG/PCA/UMAP) runs when z supplied; zero-count cell w/ z -> advice "supply z" (already supplied); tiny batch UMAP 'k capped' warning | P2 | decon.R:715-722, 957-963 | s10_umap_seed.R | skip init w/ z; keep accurate zero-count check | S/M | I-d (decision: no UMAP breaks plotDecontXContamination) |
| A1-17 | D1-11 | delta=c(0,0) accepted (cor 0.009); delta NA cryptic | P2 | decon.R:1070 | s05 | delta>0, no NA | S | I-b |
| A1-18 | D1-12 | negative counts pass through to output; .checkCountsDecon(1:10) misleading; NA check densifies DelayedArray | P2 | decon.R:867-873 | s05 | sparse-aware min>=0, anyNA, dims msg | S | I-b |
| A1-19 | D1-09 | batch length/NA unvalidated (cryptic Dimnames errors; NA -> "Missing value in 'counts'") | P2 | decon.R:443,479 | s05, s08 | as A1-01 | S | I-b |
| A1-20 | new | marker plots matrix path: single marker -> DelayedArray seed error (no drop=FALSE :214); numeric-like permuted colnames -> cluster NA | P2 | plot_decontx.R:214, 391-400 | s12 | drop=FALSE, match by name | S | I-b |
| A1-21 | new | simulateContamination ignores NRange when min==max (sample() length-1 gotcha) | P2 | decon.R:1185-1188 | s13_sim_rfi.R | avoid sample() on length-1 | S | I-b |
| A1-22 | D1-07 down | factor z unused levels -> NaN phi columns only (contam diff ~2e-16) | P3 | decon.R:500, 894 | s04, s04b | droplevels per batch | S | I-b |
| A1-23 | D1-16 | same seed per batch: equal-size batches identical init draws, but batch result == batch run alone, order-independent, RNG restored | P3 | decon.R:521-523 | s10 | recommend document + close (decision) | S | decision |
| A1-24 | D7-04 | plotDensity per-curve bandwidth (0.204 vs 0.390); original peak inflated 1.18-1.59x | P3 | plot_decontPro.R:61-67 | s16_plotdensity_bw.R | pooled bandwidth | S | I-b |
| A1-25 | new | stale decontX_b*_UMAP after re-run unbatched; groupClusters on batched labels ("b1-1") errors | P3 | decon.R:205-223 | s12 | clear old entries; doc | S | I-b/III |
| A1-26 | new | cosmetic: runParams lacks seed/iterLogLik (387-400); typos "bacth" :425, "columnsof" :1301-1302; simulateContamination(K=1) "Division by 0"; retrieveFeatureIndex(NA) matches NA row (celda_functions.R:130); by length>1 cryptic | P3 | various | s13, s17 | small fixes | S | I-b |

Not reproduced: D1-10 ddiMatrix latent only; input classes (dense int/double, dgC, dgT, DelayedMatrix(dgC), HDF5) bit-identical;
duplicated/missing dimnames, zero-count genes, non-"counts" assayName OK; decontPro dimnames kept, decont <= observed, parts sum, zero ADT row no NaN.
ExperimentSubset not installed (area 6).

Decisions: (1) column names for z/batch (II; length check now regardless); (2) dense in -> dense out?; (3) skip init w/ z -> no UMAP;
(4) keep user labels changes colData levels (celda/sctk); (5) close D1-16; (6) decontPro assay/altExp/layers (II);
(7) expose seed/iter/output_samples + error on unknown ... (ADR w/ D7-01).
Repo clean after area 1.

### Maintainer decisions (2026-10-18)
- (1) A1-14: keep sparse dgCMatrix output for dense input; delete dead coercion; document. AGREED.
- (2) A1-16/D6-04: keep computing init/UMAP when z supplied; fix docs + misleading zero-count message. AGREED.
- (3) A1-09: return z as a factor whose levels are the original labels (codes = first-appearance order for
  non-factor z, so numerics stay bit-identical; factor z keeps its level order minus unused levels).
  Pending sub-question: drop the batch prefix for user-supplied z in batched runs?
  Downstream: celda only uses decontX_clusters in an example; singleCellTK uses it as a violin grouping column
  (labels help). Must land with the A1-07 plot fix; groupClusters with user z switches from codes to labels (NEWS).
  ANY method: result$z changes from integer to factor (as.integer() recovers old values) -> NEWS.
- (4) A1-23/D1-16: document, close; no result change. AGREED.
- Env: stay on local Bioc 3.23 for Phase 0; before Phase I set up a separate devel library via an untracked
  repo-root .Renviron (R_LIBS_USER), not renv. Check .gitignore/.Rbuildignore cover .Renviron.

---

<!-- recovered from scratchpad/area2/findings.md -->

## Area 2: compiled code (agent report, 2026-10-18; D2-03 claim verified by lead: RcppExports.cpp regenerated 18:32, git clean)

Env: R 4.6.1, Rcpp 1.1.2, RcppEigen 0.3.4.0.2 (Eigen 3.4.0), StanHeaders 2.39.1, rstan 2.32.7, MCMCprecision 0.4.2.
Scripts s01-s10 + out_*.txt here; crash cases in subprocesses. Scratch copies pkgcopy/ (lib/), pkgmig/ (mig/lib/). Repo clean.

| ID | Title | Sev | Where | Repro/observed | Fix | Eff | Phase |
|---|---|---|---|---|---|---|---|
| D2-01 | _colSumByGroup_numeric no NA guard; R wrapper factor(group, levels=seq(K)) turns out-of-range labels into NA | P2 (was P1; no package caller) | matrixSums.c:50-88; R/matrixSums.R:16 | s01: num_bigNA exit 139 segfault; wrapper c(1,2,3,1,2),K=2 drops label-3 col silently + wild write (UB) | delete numeric variant (dead) | S | I-b/I-d |
| D2-02 | decontXInitialize no min(z)<1 guard | P2 (internal only) | DecontX.cpp:227-238 | s02: z=0 "Not a matrix" (header corruption); negative -> exit 139; NA -> negative length | range check like other exports | S | I-b |
| D2-03 | registration is GENERATED (package_native_routine_registration_skeleton scans R/ .Call names), not hand-edited | P3 (largely obsolete) | RcppExports.cpp:184-203 | compileAttributes on copy byte-identical; load_all re-runs it; backtick symbols identical; removing wrapper drops entry | drop ordering rule; no init.c; delete C routine with its R wrapper | S | I-d |
| D2-04 | dead compiled code: all 4 fns in matrixSumsSparse.cpp (0% cov), _colSumByGroup_numeric + wrapper, fastNormPropLog/Sqrt, nonzero | P2 | matrixSumsSparse.cpp:1-189; matrixNorm.cpp:35-117; R/matrixSums.R:4-19 | no R callers besides RcppExports.R/tests | delete; keep fastNormProp + integer .colSumByGroup | S | I-d |
| D2-05/06 | *ChangeSparse mutate caller px in place; px dims unchecked | P2 (dead) | matrixSumsSparse.cpp:109-117,165-173 | s09: mutate TRUE; oob exit 139; msg :167 says rows checks cols | deletion (D2-04) | S | I-d |
| D2-07 + A2-03 | int offsets overflow G*C>=2^31; int accumulation overflows silently, NA_INTEGER as INT_MIN; "pgroup" msg | P3 | matrixSums.c:25,34,42,74,82 | reasoned; int_overflow -> -2 vs rowsum 4294967294 | R_xlen_t, or base rowsum() in .simulateContaminatedMatrix and drop matrixSums.c | S | I-d |
| D2-08 | decontXEM 0/0 for empty cluster: NaN stays in that phi column only; zero-count cell -> NaN contam then fit_dirichlet error (R guards stop it first) | P3 | DecontX.cpp:114-118,123-127,131 | s03, s04 | droplevels per batch (A1-22); optionally skip zero-sum columns | S | I-b |
| D2-09 | calculateNativeMatrix exp/log round trip, coeffRef binary search, copy+wrap, unused var x (only compiler warning) | P3 perf | DecontX.cpp:317-339 (:321) | s05: direct form max rel diff 7.7e-15; 0.041 s vs 0.004 s at nnz 1.6M (-O2) | direct product, reuse i/p; needs documented tolerance | S | I-e (dead var I-d) |
| D2-10/D3-15 | MappedSparseMatrix x8 DecontX.cpp + x8 matrixSumsSparse.cpp; deprecated since Eigen 3.3, removed at Eigen HEAD | P2 future break | DecontX.cpp:7,76,159,207,221,248,283,327 | no warnings today; migration on pkgmig (sed to Map<SparseMatrix<double>> + compileAttributes) compiles clean, BIT-IDENTICAL (s10) | do after D2-04 deletions; no D2-03 dependency | S | I-d |
| D2-11 | RcppParallel no do-not-remove comment; ld duplicate -ltbb warning from generated rstantools Makevars | P3 | R/decontX-package.R:14 | install.log:626 | comment; leave Makevars | S | I-d |
| D2-13 | GPL-3 header in generated stanExports_shrinkage.h; no LICENSE note | P3 | stanExports_shrinkage.h:3-5 | - | LICENSE/inst/NOTICE note | S | I-d |
| D2-14/15 | delta msg "positive integers" but accepts 0/doubles/NaN; eta loop bound harmless | P3 | DecontX.cpp:48-51,124 | code read | fix msg; rest via D2-04 | S | I-d |
| A2-01 | load_all/devtools::test compile at -O0 (pkgbuild debug=TRUE): decontXEM 8.7x slower than installed; results identical | P3 process | - | s07, s10 | dev/bench must time an installed build (R CMD INSTALL to temp lib) | S | Phase 0 bench |
| A2-02 | //' roxygen in matrixNorm.cpp generates wrong help pages for unexported helpers | P3 | matrixNorm.cpp:5-10,30-35,55-60,81-85; man/fastNormProp*.Rd, nonzero.Rd | ls man | delete dead helpers; @noRd fastNormProp | S | I-d |

Clean: no RNG in hand-written C++; no threads (TBB is Stan's); Map exporter zero-copy; PROTECT balanced; NA z caught in EM; no interrupt check needed.

Perf candidates (area 5): benchmark installed builds only; EM loop ~4.4 ns/nnz, bounds-checked operator() costs ~20%
(bit-identical win); theta hoisting / ratio matrix / RcppParallel not bit-identical; calculateNativeMatrix ~10x (once per batch);
fit_dirichlet 0.1 ms/iter at 5k cells, 16 ms at 1e5 (not a hotspot); measure D1-27 per-batch copies; pursue EM only if it outweighs init.

Not reproduced: negative contamination -> NaN delta (can't happen today; pmax guard optional); all contamination 0 contrived;
upstream MCMCprecision fit_dirichlet fallback runif(length(NULL)) no-op.

Decisions: (1) D2-03 no longer a blocker; keep generated registration, symbol .Call() names (D1-25);
(2) delete dead compiled code; optionally base rowsum() replaces matrixSums.c + R/matrixSums.R (celda copies: coordinate);
(3) ~1e-14 tolerance for D2-09 rewrite in I-e; (4) legacyInit=TRUE now emits scater deprecation warnings
(librarySizeFactors, normalizeCounts -> scrapper): upstream deprecation bounds legacyInit lifetime per policy.

### Maintainer decisions (2026-10-18): all four AGREED
(1) D2-03 obsolete; generated registration stays; symbol .Call() names. (2) delete dead compiled code; replace
remaining .colSumByGroup with base rowsum(), drop matrixSums.c + R/matrixSums.R. (3) ~1e-14 tolerance for D2-09 in I-e.
(4) legacyInit: no change, let scater deprecation warnings show; ADR only when scater makes them defunct.

---

<!-- recovered from scratchpad/area3/findings.md -->

## Area 3: tests and coverage (agent report, 2026-10-18)

Artifacts here: cov.rds, tally_line.csv, zero_blocks.txt (baseline); draft/test-draft.R (~70 expectations, +2.7s);
cov2.rds, tally_line2.csv (suite + draft); after.R/after2.R; refs.rds (140 KB sample references). Repo clean.

### Coverage (covered/total lines, %)
| File | Now | After dead-code deletion | +draft measured |
|---|---|---|---|
| R/celda_functions.R | 23/76 30.3 | 30.3 | 100 |
| R/decon.R | 599/716 83.7 | 592/707 83.7 | 95.4 |
| R/decontPro.R | 68/71 95.8 | 95.8 | 100 |
| R/plot_decontPro.R | 0/117 0 | 0 | 100 |
| R/plot_decontx.R | 263/309 85.1 | 85.1 (87.7 w/ A3-12) | 97.1 |
| R/stan_helpers.R | 100 | 100 | 100 |
| src/DecontX.cpp | 134/156 85.9 | 85.9 | 100 |
| src/matrixNorm.cpp | 10/58 17.2 | 10/11 90.9 | 100 |
| R/ total | 1028/1369 75.1 | 1012/1346 75.2 | 96.6 |
| Overall covr printed | 74.21 | 77.9 | 87.2 |
| Overall excl. generated Stan (codecov denominator) | 1215/1705 71.3 | 1156/1513 76.4 | 89.9 (97.4 after del) |
Testing plot_decontPro.R + celda_functions.R alone -> R/ ~87.8%, overall ~87.6%.
Config: codecov.yml targets 80% (thresholds 1%/5%), ignores RcppExports + stanmodels + stanExports; no .covrignore,
so make coverage counts the 718-line generated .h (inflates local by ~3 pts).
Still uncovered after draft: decon.R:608-628 (DelayedArray/HDF5, A1-06), :995 (scrapper<1.5 branch), :1013-1017
(scater missing), :292-305 (covr artifact, A3-13), plot_decontx.R:424-432 (unreachable, A3-12).

### Exports
Tested: decontX SCE, ANY dense; decontPro ANY (real), SCE/Seurat (mocked); simulateContamination; plots structure.
Untested: decontX ANY dgC/DelayedArray/HDF5; batched SCE; batched+background; decontXcounts<- directly (getter dims only);
decontPro ambient paths, split layers; simulateContamination seed=NULL/scalar delta; retrieveFeatureIndex 7/8 branches;
plotDecontXContamination batch= + errors; marker percentage values; log1p/plotDots; plotDensity, plotBoxByCluster (0%).
stop()/warning(): 35 in R/, 9 hit, 6 with regexp. Untested: celda_functions.R:9,93,100,105,121w,136,143,151w;
decon.R:424,437,872,1013,1051,1063,1276w,1300; decon.R:1071 bare (:284,:288), 1175w bare (:121), 1203 bare (:122);
plot_decontx.R:36,47,173,469,475,479,494,501,511. C++: 22 guards in DecontX.cpp + matrixNorm.cpp:21 untested.

### Findings
- A3-01 P1 I-a M: no pinned reference outputs (only within-run comparisons) -> nothing to compare I-b/rowsum/I-d/I-e against.
- A3-02 P1 I-c M: targets unmet (R/ 75.1%; 26/35 stop paths untested).
- A3-03 P2 I-c S: local vs codecov denominators differ -> .covrignore (stanExports .h/.cc, R/stanmodels.R) + .Rbuildignore.
- A3-04 P2 I-c(bg in I-a) M: hollow: test-decon.R:172-188, 254-271 expect_s3_class on unbuilt ggplots; :272 unasserted;
  :190-197, 274-278 background checks dims only (no eta check); test-decontPro.R:3-20 sum identity only (no-op passes);
  :69-82 legacyInit "reproduces" asserts range only under suppressWarnings; test-decontPro.R:134 unused res; :32-59 use expect_no_error.
- A3-05 P2 I-a S: test-decon.R:85-86 sim1==sim2 (same default seed) -> batch misalignment invisible.
- A3-06 P2 I-a S: bare expectations test-decon.R:121,122,284,288; test-matrixSums.R:13,14,24,25 -> regexp before A1-10/A1-17.
- A3-07 P2 I-c S: expect_equivalent test-decon.R:112 (#55; baseline warning 1).
- A3-08 P3 I-c S: test-decontPro.R:118 CreateSeuratObject(dense) coercion warning (baseline warning 2) -> CsparseMatrix fixture;
  skip_if_not_installed("SeuratObject") redundant.
- A3-09 P3 I-c S: 5 real Stan fits on same fixture (5.9 of 11.2s); merge :3 + :64; add Stan-free .process_stan_vb_out test.
- A3-10 P2 I-a S: ambient happy path uncovered (decontPro.R:72,95,143); mock .call_stan_vb to capture dat$p.
- A3-11 P3 I-d S: decon.R:995 old scrapper runUmap API uncoverable; floor is scrapper >= 1.2.0.
- A3-12 P3 I-d S: plot_decontx.R:424-432 unreachable.
- A3-13 P3 I-d S: decontXcounts GET_FUN/SET_FUN closures (decon.R:292-305) show 0% though called -> explicit methods.
- A3-14 P3 I-c S: scater-missing stop (decon.R:1013) unmockable -> .hasScater() helper or nocov.
- A3-15 P3 I-c S: noisy output (verbose=TRUE everywhere; Stan "Chain 1:").
- A3-16 P3 I-a/I-c S: no snapshots/withr; top-level fixtures (test-decon.R:3, test-matrixSums.R:4-10); fixture mutation
  test-matrixSums.R:21; library(decontX) unneeded; snapshots skipped on Bioc (NOT_CRAN unset).
- A3-17 P3 I-a->I-d S: test-matrixSums.R tests dead code; pin simulateContamination output before rowsum swap.
- A3-18 P2 I-c S: C++ guards untested; direct calls cover all in 0.08s; add D2-02 guard.
- A3-19 P2 I-a/I-b S: no dgC/DelayedArray/HDF5 input test (A1-06 site); HDF5Array used at :620 but undeclared (D3-09).
- A3-20 P3 I-c S: maxIter=2 in 6 calls; fixture converges in 6-22 iter (~0.1s).

### I-a pinning plan (P = pin stays green; R = regression test, fails now)
A1-01/19: R z len 90, z="cluster", bad batch len/NA; P batched ref. | A1-02: P bg ref; R reversed rows identical, missing row errors,
unnamed unequal errors. | A1-03: P batched+bg+bgBatch ref + :424/:437 tests; R bgBatch mismatch, batch w/o bg errors. |
A1-04: R split layers -> mock sees all 10 cols, no warning. | A1-06: P DelayedArray(dgC)==matrix; R DelayedArray(matrix) -> DelayedMatrix. |
A1-07: P exact .calculateDecontXBarplotPercent for 1..K; R 0-based, {10,11,12}, {2,5,9}, unused levels. | A1-08: R 1-cell batch/complement, 1-col bg batch. |
A1-09: P as.integer(result$z) + contamination == ref; R char z label levels, factor keeps levels minus unused, batched "b-label", NA errors, groupClusters on labels. |
A1-10/17: regexp on :284/:288 first; tests for :1051/:1063; R each bad scalar; maxIter<iterLogLik records loglik. |
A1-11/13: R cell_type gap/NA/len; mocked Stan failure informative; unknown ... errors. | A1-12: P mocked dat$p; R reversed rows same p, missing ADT errors. |
A1-14/15: P dense -> dgC, dgC == dense; R dgRMatrix runs, dgTMatrix no warning, data.frame clear error. | A1-16: R message no longer says supply z. |
A1-18: R negative counts error. | A1-20: R single marker, permuted colnames. | A1-21: P simulateContamination ref; R NRange c(500,500). |
A1-22: R no NaN phi. | A1-24: plot structure then shared bw. | A1-25: R unbatched rerun removes old UMAPs. | A1-26: snapshots. | D2-02: R guard.

References: committed RDS in tests/testthat/fixtures/ (+ fixtures/make-refs.R at pre-change commit); 6 paths on
simulateContamination(C=100,G=60,K=3,seed=1) = 140 KB xz, identical on rerun: single+z, batched+z, bg, bg batched+bgBatch,
no z, estimateDelta=FALSE; plus SCE layout and simulateContamination default. Store contamination, z, decontXcounts, per-batch
theta/phi/eta/delta/logLik/iteration. Compare user-z paths with tolerance 1e-10 (cross-platform libm); no-z and decontPro Stan
not in strict committed comparisons. Bit-identity gate: dev/ script installing old + new commits in two temp libs, identical()
over all paths incl. no-z, legacyInit, decontPro 10x8 w/ and w/o ambient.

### Per-file plan
plot_decontPro.R 0 -> 100 (patchwork structure, layer_data, file= with withr tempfile, bandwidth);
celda_functions.R 30 -> 100 (all retrieveFeatureIndex branches, logfile); decon.R 83.7 -> 95.4 (~98 w/ A1-06, A3-11, A3-13);
plot_decontx.R 85.1 -> 97.1 (100 w/ A3-12); decontPro.R -> 100 (ambient mocked); DecontX.cpp -> 100 (direct guards);
matrixNorm.cpp -> 100. R/ ~97-98; overall excl. Stan ~97-98.

### Timing
11.2s wall for 23 blocks (load_all -O0): decontPro 5.95s, decon 5.25s, matrixSums 0.01s. Slowest: decontPro 1000-draw fit 2.54s,
legacyInit 1.22s, 3 Stan fits ~1.1s each. Draft +2.7s (1.4s legacyInit seed=NULL, droppable); merging Stan fits -1.1s -> ~12-13s.

### Decisions (agent's list)
1 .covrignore + targets on codecov denominator; 2 committed RDS refs (1e-10) + dev/ cross-commit identical() gate;
3 scrapper floor 1.5 + delete :995 or nocov; 4 delete plot_decontx.R:424-432, explicit decontXcounts methods;
5 .hasScater() or nocov; 6 HDF5Array to Suggests or drop special case; 7 Stan refresh=0 vs test-only silence;
8 regexp/class primary, snapshots for wording; 9 legacyInit test: pin stored numbers or rename.

### Maintainer decisions (2026-10-18): all nine AGREED as recommended
1 .covrignore, targets on codecov denominator. 2 committed RDS refs (1e-10, user-z paths) + dev/ cross-commit identical() gate
(no-z, legacyInit, decontPro local only). 3 scrapper >= 1.5, delete :995 branch (deps PR). 4 delete plot_decontx.R:424-432,
explicit decontXcounts methods. 5 .hasScater() helper. 6 keep HDF5-in -> HDF5-out; HDF5Array to Suggests (fold into ADR-0005).
7 silence Stan in tests only; user-facing quiet default/verbose goes in decontPro-args ADR. 8 regexp/class primary, snapshots wording.
9 rename legacyInit test; dev/ gate carries reproducibility.

---

<!-- recovered from scratchpad/area4/findings.md -->

## Area 4: code quality (agent report, 2026-10-18; A4-09 ordering claim verified by lead)

Scripts: fnlen.R, melt_identity*.R, plyr_identity.R, exact_repl.R, colsum_factor.R, lint4.R (on rcopy/). Repo clean.

### Argument drift (decontX)
Sites: SCE :140, ANY :235, .decontX :354, .decontXoneBatch :674, runParams :387, .decontxInitializeZ :902; forwarding blocks
:181-198, :269-286, :503-519, :523-539, :715-722. Drift: maxIter 500 vs 200 (.decontX, oneBatch); convergence 0.001 vs 0.01 (oneBatch);
varGenes 5000 vs NULL (core) vs 2000 (initializeZ); dbscanEps 1 vs NULL; verbose TRUE vs FALSE (.checkBackground); batch unused in oneBatch;
runParams: assayName/iterLogLik/seed missing, batch overwritten with processed vector (:447) and used as data by SCE method (:205,:211).
NOT reachable today (all internal calls pass every arg) -> D1-20 re-rated P1 -> P2 (latent).
decontPro: 3 methods same defaults (delta_sd 2e-5, background_sd 2e-6); generic dispatches on cell_type too (no signature=);
SCE assay "counts" hardcoded :70/:72; Seurat default assay :93/:95; ... swallowed; seed/iter/output_samples/pars hardcoded
(stan_helpers.R:20-23); delta_sd/background_sd documented as "Prior variance" but used as normal SD (shrinkage.stan:42-43).

### Function lengths (>50): .decontX 307 (decon.R:354-660); .decontXoneBatch 190 (:674-863); plotDecontXMarkerPercentage 126;
plotDecontXMarkerExpression 115; .simulateContaminatedMatrix 114; .processPlotDecontXMarkerInupt 97; decontX SCE 92; plotBoxByCluster 86;
retrieveFeatureIndex 80; .process_stan_vb_out 79; plotDensity 73; .decontPro 72; plotDecontXContamination 70; decontX ANY 53.

### Proposed core (all in decon.R)
- .decontXParams(maxIter, convergence, iterLogLik, delta, estimateDelta, varGenes, dbscanEps, seed, legacyInit, logfile, verbose): NO defaults;
  validates scalars (A1-10, A1-17); absorbs .checkDelta (:1069-1077), .processvarGenes (:1046-1055), .processdbscanEps (:1058-1067), logfile check
  (celda_functions.R:8-13); returns runParams (+seed, iterLogLik).
- .decontX(counts, z, batch, background, bgBatch, params): ~35-line orchestrator -> list(decontXcounts, contamination, z, batch, estimates, umap, runParams).
- .decontXPrepareInputs(...): :403-449 + .checkBackground (:1268-1315, matrices) + .checkCountsDecon once (sparse-aware); z/batch length+NA (A1-01/19),
  bg row matching (A1-02), bgBatch (A1-03).
- .decontXSubsetBatch(inp, bat): :478-501 with drop=FALSE, CsparseMatrix.
- .withSeed(seed, code): replaces :502-541, :1024-1039, :1123-1146; bit-identical.
- .decontXRunBatch -> .decontXoneBatch(counts, z, background, params) ~25 lines: .decontXInitBatch (:695-733), .decontXInitEM (:748-771),
  .decontXRunEM (:773-846; merge duplicated decontXEM calls via estimate_eta; drop stopIter/deconMethod). Keep RNG order in one seed scope.
- .decontXNativeCounts (:553-560); .decontXAssembleCounts (:409-416, :562-574; one cbind + reorder, D1-27, identical() check);
  .decontXMergeLabels (:589-594; decision 3); .decontXRestoreClass (:606-639; keep HDF5; delete dead canCoerce);
  .decontXUmaps (moved out of SCE method :205-223); .decontXLog.
Container touches: SCE method :162-171, :176, :179, :201-202, :205-223, :226, :228; ANY :255; accessors :292-347; legacy init SCE (:1019-1022, fine);
decontPro.R:70,72,93,95; plots (below). Validation in wrong layer: :689-690, :712-713, :790/:803/:908 isTRUE, :779/:816/:822, :880, :419-449.

### Findings
- A4-01 P2 M I-d: D1-20 latent drift; .decontXParams no defaults + formals-parity test across methods (more important in II).
- A4-02 P2 M I-d: validation scattered/per batch (:880 root of A1-01; delta/varGenes/dbscanEps per batch; NA check after dgC coercion :689; logfile per message;
  isTRUE coercion). I-b fixes in place first; I-d moves them.
- A4-03 P2 M I-d: D1-21 monoliths (decomposition above).
- A4-04 P2 S I-d: runParams used as data (:205, :211); stores processed batch + raw z; lacks seed/iterLogLik. Core returns batch, umap.
- A4-05 P2 S I-d/II: core logic in SCE method (UMAP assembly :205-223; bg dedup on containers :158-177; bgAssayName fallback); SCE x requires SCE bg (A1-15).
  Per-class .readCounts(obj, assayName).
- A4-06 P3 S I-d: 3 seed blocks -> .withSeed().
- A4-07 P3 S I-d: dead/redundant: duplicated decontXEM calls :781-807; inert stopIter :736-737,:780; deconMethod :693,:747; umap<-NULL :696;
  duplicate z= key :577/:585; unused batch arg :676; totalCells==nC :404/:406; 1e-20 literal x5 (:559,:759,:792,:805,:829); redundant as.factor plot_decontx.R:528.
- A4-08 P3 S I-d/lint: 17 .logMessages calls repeat logfile/append/verbose; 6 identical rule lines; cat() celda_functions.R:14 (write() changes log text: space before newline).
- A4-09 P2 S I-d deps: plyr -> match() NOT identical for character z >= 10 labels (old levels string-sorted: codes differ, phi/eta permuted, contamination 3.3e-16).
  Exact: u <- match(z, unique(z)); if (is.character(z)) u <- as.character(u); as.factor(u). Invalidates decision-3 premise "first-appearance codes".
- A4-10 P3 S I-d deps: reshape2 5 sites; faithful helper (var.convert via utils::type.convert) identical; A1-07/A1-20 fixes change plot_decontx.R sites first.
- A4-11 P3 S lint/II: plot duplication: marker assay loops :191-222 vs :369-396 (drop=FALSE only at :390); ncol defaults differ;
  plotDensity/plotBoxByCluster repeat fill scale (70-73/200-203), breaks (66-69/196-199), patchwork (94-101/224-235), ggsave (105-114/237-246); @inheritParams.
- A4-12 P3 S II: decontPro methods triplicated (decontPro.R:64-122); .decontProParams; setGeneric(signature = "filtered_counts").
- A4-13 P3 S I-d: internal naming (.decontxInitializeZ*, .processvarGenes, .processdbscanEps, ...Inupt typo, snake .call_stan_vb); mixed local var styles.
  Exported (breaking, flag only): decontPro snake_case args; plotDensity/plotBoxByCluster generic names (affy exports plotDensity).
- A4-14 P3 S lint: stan_helpers.R lacks @noRd (man/dot-call_stan_vb.Rd, dot-process_stan_vb_out.Rd); decontPro.R:51 stray #' NULL, generic no @export/@rdname;
  variance vs SD docs; vendored celda helpers no version header.
- A4-15 P3 S lint: seq(K/C/G) decon.R:1173-1252, celda :37, matrixSums :16/23 -> seq_len; sample(seq(n)) -> sample.int (:1173,:1210); sum(is.na())>0 -> anyNA
  (:868, celda :133); is(x,"list") -> is.list (plot_decontx.R:172,:493); return() celda :42,:162; dead assignments plot_decontx.R:168; embedded newline in
  error strings (:36-39, :469-470); measure.var partial matching (plot_decontPro.R:57,:179); 28 base-masking names; matrix(rep(...)) decontPro.R:176-183.
- A4-16 P3 S I-d deps: bare with_seed (:1025,:1134) vs withr:: (:521); redundant importFrom(plyr); dup useDynLib + string .Call go with matrixSums.R;
  stats (8 uses), grid (plot_decontx.R:262) used via :: but not in Imports (+utils for type.convert); rstantools template imports keep w/ comment;
  single-use imports: dbscan, plyr, reshape2, S4Vectors, patchwork, SeuratObject; DelayedArray only marker-percent rowsum/colsum (:567-573) + class restore.
- A4-17 P3 S post-I-d move PR: R/decontXcounts.R, R/simulateContamination.R (+ .rdirichlet); core stays in decon.R. (file split = structural -> ADR-0004)
- A4-18 P2 M II: plots read containers directly (below).
- A4-19 P3 S II: marker plots subset whole SCE (x[geneMarkerIndex,] :197/:375; x[,!na.ix] :524) copying all assays.
- A4-20 P3 S I-e: .process_stan_vb_out denominator x3 (stan_helpers.R:103-110), background_est assigned twice (:47,:95), unscaled_rates reused (:70,:87).

### plyr/reshape2 replacements
decon.R:895 exact replacement above (identical, 5 cases); decon.R:1182 match() identical (integer); plot_decontx.R:378, :391, :578 faithful
var.convert helper identical (incl. numeric-like/unnamed/duplicated names, layer_data); naive factor(colnames) differs = A1-20 fix;
plot_decontPro.R:57, :179 hand-built data.frame identical (layer_data); matrixSums -> unname(t(rowsum(t(cellRmat), z))) identical.

### Phase II obstacles
plot_decontx.R:30-34 metadata(x)$decontX$estimates or x$estimates (R6 $ returns NULL -> test is.list first); inherits SCE x5 (:30,:191,:369,:466,:473)
(ExperimentSubset reads top level); assay(x[idx,]) :196-199,:374-377; colnames :399; colData decontX_clusters :467-472; colData(x)[,z] :474-477;
ncol :478,:518; x[,!na.ix] :524; retrieveFeatureIndex rowData(x)[,by] (celda :98-102). Getter returns list(estimates, clusters, mats, cellMeta col,
featureSearch vector). Also runParams$batch as data; UMAP assembly in SCE method; assayName semantics per class; metadata nesting (uns spike);
decontPro generic signature, ..., hardcoded assay; naming mix for new decontPro args.

### Decisions (agent's list)
1 params pattern A (full signatures + no-default constructor + formals-parity test) vs B (...); 2 plyr exact (internal order) vs first-appearance;
3 reshape2 sequencing (I-b drops melt in plot_decontx.R first); 4 core return shape (+batch, +umap, runParams = validated params; raw z?;
ANY returns extra fields); 5 pure-move PR; 6 cat() fix in decontX vs celda first; 7 naming of new decontPro args, renames on clean break (default none);
8 decontPro signature="filtered_counts"; declare stats/grid/utils in Imports.

### Maintainer decisions (2026-10-18)
1 pattern A YES. 2 plyr exact replacement (internal order bit-identical) YES. 3 reshape2 sequencing YES. 4 core return shape YES
(+batch, +umap; runParams = validated params +seed +iterLogLik; keep raw z; ANY returns extra fields, NEWS). 5 pure-move PR YES (record in ADR-0004).
6 fix cat() in decontX + celda issue YES.
7 decontPro naming: BOTH. camelCase are the primary arguments; existing snake_case names stay as aliases that map onto them
  (design in ADR-0005: likely permanent aliases, error if both given; dispatch-arg alias via a "missing" method).
8 asked what it does (explained: generic currently dispatches on filtered_counts AND cell_type; restrict to filtered_counts;
  needed for the alias design; Imports declarations are no-behaviour-change) -> pending confirmation.
Label ordering (decision 3 revision): internal codes keep today's exact order (bit-identical); returned factor levels in natural
order (user factor levels, else factor(z) order); phi/eta columns named by cluster label; as.integer(result$z) may change -> NEWS.
Taken as YES under "all others are yes" (confirm in summary).

---

<!-- recovered from scratchpad/area5/findings.md -->

## Area 5: speed and memory (agent report, 2026-10-18; dense-output claim corrected by lead)

Script: dev/bench/bench-decontx.R (Rscript dev/bench/bench-decontx.R [sim|real|profile] [outdir]); installs a copy of the working tree
via R CMD INSTALL into <outdir>/lib (never load_all); each case in a fresh Rscript worker under /usr/bin/time -l for peak RSS; env vars
DECONTX_BENCH_LIB/LABEL/SIZES/RUNS/PRO_CELLS/ITER/PROFILE_CELLS/PBMC; lint clean. Outputs here: out-sim/, out-real/, results-10k-50k-pro5000.csv,
targeted*.txt, profreal.R. Installed -O2 build throughout. Sim: G=2000, K=10, NRange 1000-4000, seed 12345, dgC; 4 contiguous batches;
background = as many empty droplets (10-200 UMIs); decontPro 40 ADTs x 5 types.

### Results (RSS incl. ~0.9 GB package-load baseline)
| cells | run | wall s | R alloc MB | peak RSS MB (before) | EM iters |
|10k|z_supplied|7.56|261|1669 (967)|11|
|10k|z_estimated|7.60|261|1516 (961)|11|
|10k|batched4|7.40|644|1884 (1043)|12;12;12;13|
|10k|background|7.62|264|1677 (1054)|13|
|50k|z_supplied|26.3|1236|2913 (1372)|12|
|50k|z_estimated|27.2|1249|2910 (1375)|14|
|50k|batched4|36.0|3087|2986 (1311)|12;11;11;11|
|50k|background|26.7|1250|2888 (1290)|14|
|100k|z_supplied|49.0|2415|3711 (1678)|11|
|100k|z_estimated|53.1|2471|3929 (1666)|16|
|100k|batched4|59.3|6063|4559 (1674)|11x4|
|100k|background|50.2|2441|3829 (1660)|13|
|pbmc4k 4340x19773 nnz 5.7M|dgC z est|5.84|463|1563 (1061)|90|
|pbmc4k|SCE HDF5|18.7|2295|2351 (1356)|90|
|decontPro 40x2000||398|n/a|1468 (915)|12,700 ADVI|
|decontPro 40x5000||260|n/a|1584 (959)|2,500 ADVI|
Simulations converge in 11-16 EM iterations vs 90 on PBMC: simulations understate EM; use real mode for EM PRs.

### Phase breakdown
10k z_sup 7.6s: UMAP 6.7 (89%), PCA 0.76, EM ~0.3, native 0.18. 50k z_sup 26.1s: UMAP optimize 15.2 (58%), kNN 3.7, PCA 3.9, EM 1.3 (5%),
native 0.93, assembly/restore 0.43, loglik 0.15, writeback 0.04, fit_dirichlet ~0. 50k batched4 35.8s: UMAP 25.7 (per-batch epochs), PCA 3.7,
kNN 2.9, assembly 1.1, EM 0.9. 100k z_sup 49s: UMAP 36.7 (75%), PCA 8.4, EM 0.22 s/iter, native 1.85, assembly 0.41 (1.48 batched), dbscan 2.7.
pbmc4k dgC 5.4s: UMAP 2.5 (47%), EM 1.8 (33%; 90 iter, 19 ms/iter, 3.4 ns/nnz), PCA 0.43. pbmc4k HDF5 17.2s: as(,"HDF5Matrix") restore 10.8 (63%).
Projection 100k real-like: init ~45s, EM ~21s, native 1.9s.
Ranked: 1 UMAP+kNN (runs with user z), 2 EM on real data, 3 PCA, 4 HDF5 restore, 5 calculateNativeMatrix, 6 batched assembly, 7 dbscan.

### Memory
Single batch 100k full copies: reorder estRmat[,allCol] 652 MB; native result 652 MB + transient C++ copy; cbind 652 MB; is.na(counts) 217 MB.
Batched4: 22 subscript allocs (3.26 GB) + 8 cbind (1.63 GB) ~7.5 copies; single cbind + one reorder identical, 0.64 vs 1.48 s, ~2 copies.
Peak RSS delta ~3x dgC size, peak in assembly. library(decontX) ~0.9 GB RSS, 2.9 s (SCE stack 470 MB, rstan 170, SeuratObject 110).
Dense input -> dgCMatrix output (lead verified on installed build; agent's dense-restore estimate was wrong).
HDF5 input: result written to session HDF5 dump dir (tempdir) -> saved SCE points at a temp file that disappears.

### Candidates
- A5-01 init/UMAP with user z: ~90% runtime at 50-100k w/ z, ~55% pbmc4k; bit-identical numbers, output change; S, risk M; I-e + API. DECISION.
- A5-02 opt-in numThreads (default 1): UMAP 4 thr 36.7->24.1 s (100k), bit-identical; PCA 8.4->3.0 s but ~3e-11 diff; S-M; I-e. DECISION.
- A5-03 unchecked indexing in EM loop: ~20% of loop; ~4 s at 100k real-like; bit-identical; S; I-e. Only EM change worth doing.
- A5-04 calculateNativeMatrix direct form + reuse i/p: ~1.7 s at 100k, -0.65 GB peak; 1e-14 tolerance (agreed); S; I-e.
- A5-05 D1-27 assembly: single fast path + batched one cbind/reorder: 0.4 s & ~1.3 GB peak (single), 0.9 s & 3.3 GB alloc (batched); identical; S; I-d.
- A5-06 HDF5 restore: 63% runtime on pbmc4k HDF5 + temp-file-backed result; values identical, class changes; S, risk M; I-d/II. DECISION.
- A5-07 anyNA(counts@x): 217 MB, 0.1-0.2 s at 100k; identical; S; I-b/I-d.
- A5-08 decontPro sample_file + unlink: TMPDIR CSV ~90 B per ADTxcell (1.8 GB at 200x100k); identical; S; I-e.
- A5-09 SeuratObject to Suggests: ~110 MB at load; M; II; low priority.
Not pursued: fit_dirichlet, loglik frequency, dbscan, RcppParallel EM, ratio matrices.

### decontPro
98% Stan ADVI gradient; rstan post-processing ~2%; prep/process ~0.01 s. Gradient 16 ms (2k cells) / 44 ms (5k): O(N*M)/iter.
ADVI iterations unpredictable (12,700 at 2k, 2,500 at 5k; cap 50,000 -> worst ~80 min at 5k). Peak RSS +0.55-0.62 GB: #46 fixed at this scale.
No bit-identical speedups; real gains need VB settings/model changes (non-goal or ADR).

### Maintainer decisions (2026-10-19): all AGREED
1 A5-01 REOPENS area-1 decision 2: keep init/UMAP by default (unchanged output) + new opt-out argument to skip init when z supplied;
  plotDecontXContamination accepts any existing reducedDim (e.g. user "UMAP").
2 A5-06 REVERSES area-3 decision 6: HDF5 input -> return in-memory sparse result wrapped as DelayedArray (class stays DelayedMatrix,
  values identical); users save to disk explicitly (saveHDF5SummarizedExperiment); HDF5Array NOT added to Suggests; NEWS.
3 opt-in numThreads (default 1); document ~1e-11 PCA differences for numThreads > 1.
4 decontPro: no model / VB-default changes; iter exposed via decontPro-args ADR; document runtime and cap.
5 real (PBMC) mode is the reference for EM-focused I-e PRs.

---

<!-- recovered from scratchpad/area6/findings.md -->

## Area 6: object coupling (agent report, 2026-10-19; AnnData dispatch verified by lead)

Scratch lib area6/lib: ExperimentSubset 1.22.0, TreeSummarizedExperiment 2.20.0, anndataR 1.2.2 (Bioc 3.23). Global: rhdf5 2.56.1,
SpatialExperiment 1.22.0, SeuratObject 5.4.0, zellkonverter 1.22.0. Scripts s01-s14, verify_dispatch*.R. Repo clean.

### ExperimentSubset
Today: decontX/decontXcounts/decontPro on Subset{SCE,Spatial,Tree}SE dispatch to SCE methods; decontX(es) runs, all results top level,
equal to matrix method, subsets survive. subsetName= -> unused argument. assayName="sub1" -> ES assay() resolves subset, core runs on 60x70,
then colData write decon.R:201 errors "60 elements in value to replace 100" (A6-01). decontPro(es) uses top-level counts. decontXcounts(es)
before a run returns NULL (ES .assay() returns NULL for unknown names once subsets exist).
Spike: assay(es, "sub1") = parent sliced, dimnames/class kept; works via SummarizedExperiment::assay without attaching ES; is() works unattached.
Names NOT unique: lookup top-level -> subset names -> internal assays, last subset wins; a subset assay "decontXcounts" is returned by top-level
decontXcounts(es). subsetColData(es,s)<-cd works, keeps factors; ES bug: getter with parentColData=FALSE returns rownames 1..n and writing back
renames cells; setter accepts mismatched rownames -> set rownames(cd) <- colnames(assay(es, s)). reducedDim(es, name, subsetName=) get/set;
listing needs ExperimentSubset::reducedDimNames (masks SCE's when attached). metadata(es, subsetName=) get/set works ("nope subset does not exist"
for unknown). altExp(es, e, subsetName=) exists. es[i,] leaves stale subset indices (A6-07).
Maintenance: 1.22.0, maintainer Irzam Sarfraz, Bioc since 3.12; last NEWS 1.1.0 (2021-05); only version bumps since; no URL/BugReports.

### Seurat (SeuratObject 5.4.0)
Today: decontX(seu) -> ANY -> "no method ... coercing 'Seurat' to 'dgCMatrix'"; decontXcounts(seu) no method; decontPro split-layer bug confirmed
(4 of 8 cells, warning only); plotDecontXContamination(seu) errors.
Prototype works: Assay5 JoinLayers on extracted copy (user object untouched), LayerData counts columns == Cells(seu) (interleaved batches checked);
v3 Assay has no JoinLayers -> branch on inherits(a, "Assay5"). Split layers -> batch: Layers(a, search="counts") + sub("^counts\\.","") + Cells(a, layer=)
reproduces the split exactly (dotted/underscored names OK). seu[["decontXcounts"]] <- CreateAssay5Object(counts=dgC) works (non-integer kept,
permuted cells realigned, subset of cells accepted, renamed cells error, features need not match RNA). meta.data: factors kept; name vectors by Cells.
Reductions: key must be alnum + "_" -> key "decontXUMAP_", name "decontX_UMAP"; per-batch keys sanitized + collision check; NA-padded or batch-only OK.
Misc(seu,"decontX") round-trips nested list with factors; re-run warns "Overwriting miscellanous data". validObject TRUE.

### AnnData (anndataR 1.2.2)
Today: decontX(adata) errors (coercing AnnDataView); decontPro errors.
Dispatch: setOldClass("AbstractAnnData") alone does NOT dispatch (falls to ANY) [lead verified]; full vectors e.g.
setOldClass(c("InMemoryAnnData","AbstractAnnData","R6")) dispatch [lead verified]. Classes: InMemoryAnnData, HDF5AnnData, ZarrAnnData (1.1.3+),
ReticulateAnnData, AnnDataView; none exported (A6-03). Note: registering "AbstractAnnData" alone first blocks later full-vector registration.
X/layers cells x genes; t(X) genes x cells dgC with dimnames; layers accept dgC/dgR/dense/DelayedMatrix; wrong orientation errors clearly;
obs factors fine; obsm NA rows fine. In place: ad2 <- f(ad) modifies ad. AnnDataView refuses writes (hint as_InMemoryAnnData()); backed HDF5AnnData
mode "r" write fails cryptically ("HDF5. Links. Can't delete message."); mode "r+" writes go to disk immediately.
h5ad round trip of uns: nested lists, factors, NULL, logicals, scalars, matrix values survive; length>1 vectors -> unnamed 1-d arrays; matrix dimnames
dropped; list order alphabetized. Layers, obs factors, obsm round-trip. Consequence: plotDecontXContamination(rt$uns$decontX) fails at render
(.data[[NULL]] from plot_decontx.R:67-68) (A6-04). h5ad IO needs rhdf5 (anndataR Suggests); no Python.

### SCE / altExp / SPE
assayName="raw" works; no "counts" -> clear error. decontPro(sce, ct, altExpName="ADT") today swallows the arg and decontaminates RNA (A4-12).
altExp prototype (read assay(altExp(sce,"ADT"),"counts"), .decontPro, write assay "decontPro" + metadata(ae)$decontPro, altExp(sce,"ADT") <- ae) works.
decontX(spe) dispatches to SCE method, returns SPE with spatialCoords/sample_id intact; plots work. Batched SCE writes decontX_<batch>_UMAP (decon.R:217).

### Corrected read/write table
| Class | Input | Counts out | Per-cell | UMAP | Run info |
| SCE/SPE | assay(x, assayName) validated vs assayNames; decontPro assay(altExp(x, altExpName), assayName) | assay decontXcounts; decontPro assay "decontPro" in altExp or main | colData | decontX_UMAP or decontX_<b>_UMAP (NA-padded) | metadata(x)$decontX; decontPro in altExp/main metadata |
| ES Subset{SCE,Spatial,Tree} | SummarizedExperiment::assay(es, subsetName) | setSubsetAssay(es, s, m, unique name e.g. "<s>_decontXcounts") | subsetColData(es,s)<-cd w/ rownames = subset colnames | reducedDim(es,"decontX_UMAP",subsetName=s)<- | metadata(es, subsetName=s)$decontX |
| Seurat | Assay5: JoinLayers on copy + LayerData counts; v3: LayerData; Cells order | seu[["decontXcounts"]] <- CreateAssay5Object(counts=) | meta.data named by Cells | reduction decontX_UMAP key decontXUMAP_; per batch sanitized | Misc(seu,"decontX") cleared first |
| AnnData concrete (not views) | t(X) or t(layers[[l]]), dimnames obs_names/var_names | layers[["decontXcounts"]] <- t(m) | obs | obsm[["X_decontX_umap"]], per batch X_decontX_<b>_umap | uns$decontX (lossy, A6-04) |
Getter .decontXResults(x, subsetName): estimates (per batch contamination + n x 2 UMAP, index columns by position), clusters, cellMeta,
counts(assayName) mapping logical "counts"/"decontXcounts", featureMeta (rowData / meta.features / var), reducedDim(name); plots extract matrices then subset rows.

### Dependencies (local evidence only, no web)
anndataR 1.2.2 (2026-09-18), Bioc since 3.22 (Oct 2025), maintainer Robrecht Cannoodt (scverse), Imports 9, 25 recursive, 0 new beyond decontX's tree,
patch backports 1.2.1/1.2.2, BugReports GitHub; h5ad needs rhdf5. ExperimentSubset 1.22.0, since 3.12, Depends SE/SCE/SPE/TreeSE, 87 recursive, 17 new
(magick, treeio, ape, dbplyr, RSQLite, BiocFileCache, Biostrings, ...), last NEWS 2021. zellkonverter 1.22.0, maintainer Luke Zappia, Imports 11 incl.
basilisk/reticulate, 41 recursive, 3 new + conda Python env at runtime, returns SCE (no in-place AnnData). -> Suggests for ES; anndataR over zellkonverter.

### Findings
- A6-01 P3 I-b S: ES assayName=subset name -> core runs on subset, colData write errors (decon.R:179, :201); typo assayName -> NULL -> cryptic. Validate assayName %in% assayNames(x), hint subsetName.
- A6-02 P2 II S: is(x,"SubsetSingleCellExperiment") FALSE for SubsetSpatial/SubsetTree -> .isExperimentSubset() via .hasSlot(x,"subsets") or all three classes.
- A6-03 P2 II S: setOldClass("AbstractAnnData") does not dispatch; register 5 full class vectors + inherits(x,"AbstractAnnData") redirect in ANY methods.
- A6-04 P2 II (plot fix I-b w/ A1-07) S: contamination plot uses UMAP colnames (plot_decontx.R:59, :67-68) -> fails after h5ad round trip or unnamed reducedDim; uns drops names/dimnames. Name columns by position; AnnData getter rebuilds estimates from obs/obsm/names; parity tests compare as.vector after round trip.
- A6-05 P2 II S: ES subset-assay names not unique, last wins, top-level decontXcounts(es) can return a subset's matrix, missing -> NULL. Unique generated name + collision check; decontXcounts(es, subsetName=) errors on NULL; ES issue.
- A6-06 P3 II S: ES subsetColData rownames bug -> set rownames before write; ES issue.
- A6-07 P3 II S: es[i,] stale subset indices -> getter-based plots extract matrices first; ES issue.
- A6-08 P3 II S: Seurat notes: v3 no JoinLayers; reduction key sanitizing/collisions; Misc<- warns on rerun; unnamed meta.data positional.
- A6-09 P2 I-b/II: confirms decontPro Seurat split-layer bug (A1-04) and swallowed altExpName (A4-12).
- A6-10 P3 II S: AnnDataView refuses writes; backed mode "r" cryptic -> clear errors (as_InMemoryAnnData(), mode="r+"), warn r+ writes to disk.
- A6-11 P3 II S: per-batch UMAP names missing from plan table -> specify per class in ADR-0005.

### Decisions (agent's list)
1 ES subset assay name "<subsetName>_decontXcounts" + decontXcounts(es, subsetName=); 2 ES run info metadata(es, subsetName=s)$decontX;
3 decontPro on ES: subsetName on main assay only, defer altExp-in-subset; 4 AnnData: 5 setOldClass + inherits redirect; 5 AnnDataView -> error w/ hint;
backed "r" -> error w/ mode="r+" hint; 6 uns: run params + per-batch scalars/phi/eta only, per-cell in obs/obsm, getter rebuilds; 7 Seurat separate
decontXcounts assay; per-batch reductions batch-only vs NA-padded; 8 plot assayName logical names mapped by getter; 9 decontPro SCE altExp output:
assay decontPro in altExp + metadata(altExp); 10 rhdf5 in h5ad test skips; ES Suggests-only; raise A6-05/06/07 as ES issues.

### Maintainer decisions (2026-10-19/20)
- Area-6 decisions 1-3, 7-9 AGREED as recommended (ES: subsetName accessor, run info in subset metadata, decontPro subsetName main assay only;
  Seurat separate decontXcounts assay, batch-only per-batch reductions; plot assayName logical names; decontPro altExp output in altExp).
- ExperimentSubset: fix upstream FIRST, no workaround code in decontX. Filed campbio/ExperimentSubset #2 (assay name scoping + NULL for unknown),
  #3 (subsetColData rownames), #4 (es[i,] stale indices), #5 (common virtual class/predicate). Repro: area6/es_issues_repro.R, bodies in area6/issues/.
  Consequences: write plain "decontXcounts" inside each subset once #2 scopes names (no generated names/collision checks); is(x,"ExperimentSubset")
  once #5 lands (no helper); no rownames reset once #3 lands. Keep: assayName validation (A6-01, I-b); plots extract matrices first (A4-19).
  ES support moves to END of Phase II, after ES devel fixes; Suggests floor on fixed ES version.
- AnnData: thin anndataR wrapper (Suggests), same class in/out, in place. NOT via a temp SCE: extract the matrix (t(X) or t(layers[[l]])),
  resolve z/batch (vector or obs column), call the container-agnostic core .decontX() (which runs the quick clustering when z is NULL),
  write results back to layers/obs/obsm/uns -- the same per-class pattern as Seurat/ES (plan section 5). Dispatch via inherits(x,"AbstractAnnData")
  redirect in ANY methods (no setOldClass). Views: let anndataR's clear error through. Read-only backed: one check with mode="r+" hint if the mode is
  queryable, else document. uns: run params + per-batch estimates; per-cell in obs/obsm. Plot getter reads obs/obsm/layers. Replaces area-6
  decisions 4-6; decision 10 (rhdf5 skip) stays. zellkonverter rejected (no AnnData object in R; Python via basilisk for writing).
  Verified the h5ad->SCE->h5ad route via anndataR works (s15_sce_roundtrip.R) -> document as an alternative.

---

<!-- recovered from scratchpad/area7/findings.md -->

## Area 7: docs (agent report, 2026-10-20)

Renders on O2 install of a git ls-files copy (area7/lib), R 4.6.1, Bioc 3.23, macOS arm64. Scripts: render.R, time_ex.R, proto_pro.R,
pbmc_time.R, seurat_check.R; knitted md/*.md. Repo clean.

### decontX.Rmd
Copy-verbatim breaks: L163-177 cellTypeMappings built in hidden echo=FALSE chunk (lapply+tapply+closure), shown code uses it; L160 falsely says
decontX figures out the mapping. L283 SingleCellExperiment() not attached (decontX only Imports SCE). L301 GetAssayData(slot=) defunct in
SeuratObject 5.4.0. L301/L304 seuratObject vs seuratObj. L298 "must be rounded" (Seurat accepts non-integer; soften).
Correctness: L118 + decon.R:105 "percentage" but proportions (0-1). L215 "groupCluster" (param groupClusters) and hardcoded c(2,6) NK+T vs actual
mapping T=2,B=4,Mono=1,NK=5 (L181/184 use 5/6). L259 $estimates$all_cells$delta undocumented (D6-07). L102 background: never says gene rows must
match (#43). L94 legacyInit note OK but long -> FAQ. No batch, doublets, citations.
Style: nesting/lapply closures L134-135, L139-143, L238-240; masking-order plumbing comment L128-130; non-default args L113-114, L265 (justified);
redundant eval=TRUE L74, L88, L237, L264; unlabeled final chunk L310; 38 lines >80; 7/25 chunks eval=FALSE; typos L72, L118, L147, L204, L253.
Unexplained: metadata(sce)$decontX, decontX_UMAP name, non-integer counts. library(scater) L112 leaks startup messages.

### decontPro.Rmd
L14 IndexEntry "decontPro"; L29 grammar; L56-62 CITEseq("cord_blood") downloads 23.3 MB, only 0.2 MB used (modes=c("scADT_Counts",
"coldata_scRNAseq") fetches 0.3 MB; "scADT_Counts" alone errors); `counts` reused for 3 objects (L61/L107/L115); L68 dim()[2]; L78-85 dplyr +
5-step pipe; npcs=10/dims=1:10 on 13 ADTs prints 6 warnings; coercion message; L89-97 RunUMAP warning, FeaturePlot unexplained; L99
as.integer(Idents()) silent 0->1 based conversion (D6-13); L106-110 duplicate eval=FALSE chunk; as.matrix on a matrix; L112 priors direction
not explained; docs say "Prior variance" (decontPro.R:12,14) but they are sd/scale of half-normal hyperpriors on tau_a/tau_b (shrinkage.stan);
L114-120 ~45 lines Stan output in HTML ("EXPERIMENTAL ALGORITHM", "MAY BE DIVERGING"); L122 only decontaminated_counts explained (ambient_counts,
background_counts, parameters (8 elements) never shown; r_est KxN, delta_est/background_est NxM lack dimnames). Missing: ambient_counts workflow,
runtime/memory, reproducibility (seed 12345, iter 50000), FAQ, citation. 11/13 chunks unlabeled; 6 lines >80.

### Build time (local; builders likely 2-4x slower; everything cached)
decontX 27.3 s: load_10X 4.7, decontX 6.6 (EM #1), newDecontX 6.8 (EM #2, init re-runs), plot_feature 1.4, library 3.0. TENxPBMCData pbmc4k EH1611-1613 13.9 MB.
decontPro 26.9 s: CITEseq 5.2, Seurat pipe 1.3, UMAP+plots 3.0, decontPro fit 12.6 (1000 cells x 13 ADTs), plots 1.3. cord_blood 23.3 MB (0.3 needed).
PBMC4k decontX 6.3 s (90 iter); with z 6.4 s (init still runs); 2000 cells 2.8 s; 1000 cells 1.5 s.
decontPro default/larger priors: 300 cells 5.4/3.0 s; 500 7.7/5.5; 1000 21.0/11.9. Real-data invariants: decont <= counts, parts sum to 4e-12.
cord_blood has NO empty droplets (all 7858 called cells, min 275 ADT counts). colData(dat)$celltype: 8 types, 7700 non-NA -> no Seurat/dplyr needed.

### decontPro rewrite outline
1 model (native/ambient/background, cite Yin 2024) | 2 load with modes= (0.3 MB) | 3 cell types: bundled celltype, drop NA, set.seed + sample 500,
as.integer(factor(labels)) explained (1..K), granularity advice (#40), own clustering pointer eval=FALSE or article | 4 run defaults, results='hide'
until quiet arg | 5 outputs: names, dims, parameters, sum self-check | 6 tuning: second fit larger deltaSd/backgroundSd, plotDensity before/after,
direction | 7 ambient_counts: only row sums used, simulated empty droplets (rmultinom), template for raw 10x; real data in III-c article |
8 runtime/memory/reproducibility table | 9 FAQ (#40 corrected D7-03 + bandwidth, 0-based Idents, priors, large data #46, granularity) | 10 sessionInfo.
~20-25 s local, 3 fits. III-b: container return, camelCase+aliases, altExpName / Seurat ADT assay.

### decontX vignette keep/cut/add
Keep: intro, install, TENxPBMCData (explain dgC), one decontX run, contamination UMAP, marker barplots, violin, cluster/prior prose, sessionInfo.
Cut/change: second EM (z + init opt-out, or eval=FALSE); hidden mapping -> shown hand-built cellTypeMappings (or labels after A1-09);
collapse L134-143 normalization to plain plotReducedDim calls; legacyInit note -> FAQ; fix Seurat section now (LayerData, library(SCE), one name),
native method in III-b.
FAQ: #40 (D7-03), #43 background rows (background[rownames(x), ]), #44 decontX before doublet detection (verify vs Yang 2020), #45 batch + z
labels "batch-label"; z factor; HDF5 -> in-memory DelayedArray; numThreads; init opt-out. Cite Yang 2020.
Data: keep real PBMC4k (~20 s after cut); if must shrink, subsample 2000 cells; no simulated data; no precompute.

### README/citation
README.md in sync (D6-09 resolved). Lacks citation section, usage snippet, Suggests note for build_vignettes; uses require().
DecontPro year: doi:10.1093/nar/gkad1032, NAR 52(1):e4, online 2023-11-17, issue 2024-01-11 (Crossref) -> 2024 (README, DESCRIPTION, CITATION agree);
AGENTS.md:23 "2023" wrong. inst/CITATION cites both correctly. decontX-package.Rd cites only RStan 2.21.7 (stale). Vignettes cite neither.
DESCRIPTION Description lacks <doi:>, starts "This package contains", Title not title case.

### Examples (local s)
decontPro 3.44 (100x10 random, no set.seed, positional priors; 10x8 = 0.95 s); decontX 2.23 (decontX 0.24, rest scater+ggplot);
plotDecontX* 1.67-1.73 each (shared example); plotDensity/plotBoxByCluster 0.17; retrieveFeatureIndex 0.00 (warns); simulateContamination 0.00;
none: decontXcounts, decontX-package; internal pages: dot-call_stan_vb, dot-process_stan_vb_out, fastNormProp*, nonzero.
@param gaps: decontX assayName, estimateDelta defaults; retrieveFeatureIndex by/exactMatch, link{regex} (celda_functions.R:74); decontPro
"Prior variance", ambient "NxM" (row sums only), cell_type 1..K rule; plotDensity/plotBoxByCluster file gets ".pdf" appended, no Default NULL,
wrong return; simulateContamination delta doc c(1,5) vs code c(1,10) (decon.R:1097-1101 vs :1120), @return "nativeMatirx"/"observedMatrix"
vs nativeCounts/observedCounts, omits z/eta/phi/markers/contamination/NByC/numMarkers. decontX @return estimates undocumented; z duplicated in
estimates$all_cells (= A4-07 duplicate z= key). decontPro @return vague. "Used only when z is not provided" false (A1-16/A5-01).

### Findings
A7-01 P1 III-a S decontX.Rmd:163-177 hidden mapping | A7-02 P1 III-a S decontX.Rmd:283,301,304 Seurat section broken | A7-03 P1 III-a S
decontPro.Rmd:114-122 Stan noise + 1/4 outputs explained | A7-04 P1 III-a M decontPro gaps -> rewrite | A7-05 P2 III-a S decontPro.Rmd:56-100
Seurat/dplyr/pipe/warnings/23 MB | A7-06 P2 III-a S decontX.Rmd:264-271 second EM | A7-07 P2 III-a M no FAQ | A7-08 P2 III-a S percentage vs
proportion (decontX.Rmd:118, decon.R:105) | A7-09 P2 III-a S nesting/closures | A7-10 P2 III-a S simulateContamination roxygen (decon.R:1097-1110)
| A7-11 P2 III-a S decontPro roxygen (decontPro.R:5-21) | A7-12 P2 III-a S no paper citations in vignettes/package Rd | A7-13 P2 I-b S
"Used only when z is not provided" (decon.R:67-90) | A7-14 P3 S AGENTS.md:23 2023 -> 2024 | A7-15 P3 III-a S decontXcounts example |
A7-16 P3 I-d S internal man pages | A7-17 P3 III-a/I-b S plotDensity/plotBoxByCluster file ".pdf" | A7-18 P3 III-a S retrieveFeatureIndex roxygen |
A7-19 P3 III-a S chunk labels, long lines, typos, IndexEntry | A7-20 P3 III-a S README | A7-21 P3 III-a S DESCRIPTION Title/Description |
A7-22 P3 I-d S z duplicated in estimates | A7-23 P3 III-a S examples (decontPro example -> 10x8).

### Decisions (agent's list)
1 PBMC4k full vs 2000 subsample; 2 second EM: wait for init opt-out vs eval=FALSE now; 3 decontPro bundled celltype, drop Seurat clustering;
dplyr out of Suggests?; 4 simulated empty droplets in vignette, real in III-c article; 5 decontPro subsample 500 vs 300; 6 Stan noise
results='hide' now vs wait for quiet arg; 7 fix decontX Seurat section now as stopgap; 8 citation 2024 everywhere + AGENTS.md fix.

### Maintainer notes (2026-10-20)
- Issue #42 (Xenium/spatial): OUT OF SCOPE for this program; maintainer plans a new spatial model after the software improvements.
- The 8 area-7 recommendations: not yet explicitly answered (ask again in the area 8 summary).

---

<!-- recovered from scratchpad/area8/findings.md -->

## Area 8: targeted methods check (agent report, 2026-10-20). Scripts s01-s05 here. #42 skipped (out of scope). Nothing posted.

### #44 doublets -> FAQ; REFUTES 2026-09 draft order (A8-03)
Model: phi_k = normalized expected native counts of cluster k; eta_k = normalized pool of other clusters' native counts (paper: eta_k = sum_{k'!=k} w_k' phi_k');
theta_j ~ Beta(delta), delta via fit_dirichlet; with background, eta fixed empirical.
Paper (Yang 2020, PMC7059395): Scrublet doublets had higher contamination (p<2e-16); all cells >70% contamination were Scrublet doublets;
contamination "orthogonal information for doublet detection"; no ordering guidance; both run on raw data (as singleCellTK runCellQC does).
Simulation (900 singlets + 100 heterotypic doublets, K=4, G=500, 2 seeds):
- doublets clustered with a parent: doublet contamination median 0.53; second parent's markers kept 0-1% (doublet looks like singlet);
  singlet MAE 0.002 -> 0.005.
- doublets in own clusters: contamination 0.09-0.10, both parents' markers kept 95%; singlet bias +0.004.
Detector AUC raw vs decontX counts (parent regime): scDblFinder 1.00 -> 0.60 (3/100 called); bcds 0.97-0.98 -> 0.76-0.80; cxds 0.98 -> 0.98
(binarizes >0; tiny non-zeros kept) but rounded 0.43. Own-cluster regime: scDblFinder on decontX 0.99-1.00. decontX contamination as doublet
score: AUC 0.998 (parent regime), 0.61 (own clusters).
DRAFT FAQ: "Should I remove doublets before running decontX? Run doublet detection (for example scDblFinder or scds) on the original filtered
counts, not on decontaminated counts, and remove the doublets before running decontX if you can. DecontX treats each cell as a mix of its own
cluster's profile plus contamination from the other clusters. When a doublet is grouped with one of its two cell types, decontX labels the
other cell's transcripts as contamination: the doublet gets a high contamination estimate (often around 50%), and its decontaminated counts
lose the second cell's markers, so a doublet detector run afterwards will likely miss it. Doublet detectors are designed for raw integer
counts; decontaminated counts are non-integer and keep many tiny non-zero values. If you already ran decontX on all cells, run the doublet
detector on the original counts. Treat very high decontX contamination as extra evidence of a doublet: in Yang et al. (2020), every cell
with more than 70% estimated contamination was also called a doublet by Scrublet."

### #40 -> reply + close; FAQ III-a; bandwidth fix I-b; invariant tests I-a
Corrected D7-03 holds: decontPro parts = counts x rate / sum of three rates, rates >= 0 (stan_helpers.R:103-110, shrinkage.stan constraints);
decontX multiplies counts by normp in [0,1] (DecontX.cpp:338-339); max(decont - counts) = 0 single/batched/background (s03); decontPro parts sum
4e-12 on real data; after CPM own-cluster markers rise in 100% of entries. LLM explanation in the issue image (prior added to native counts) wrong.
DRAFT FAQ: "Why do decontaminated values look higher than the original values? No decontaminated count is ever larger than the original
count. decontX multiplies each count by the estimated probability that it is native (a number between 0 and 1). decontPro splits each count
into native, ambient and background parts in proportion to their estimated rates, and decontaminated_counts + ambient_counts +
background_counts equals the input counts exactly. Two things can make results look higher. (1) After normalization (library size, log, or
CLR), a retained marker can rise because removing contamination shrank the cell's total; this is the intended gain in signal-to-noise, not
added counts. (2) In plotDensity(), each curve currently gets its own smoothing bandwidth. Decontamination squeezes the negative population
toward zero, so the decontaminated curve is smoothed less and its peaks look taller and spikier even though no value increased (a fix using
one shared bandwidth is planned). To check: max(decontXcounts(sce) - counts(sce)) # 0 or less; max(out$decontaminated_counts - adt_counts)
# 0 or less; all.equal(out$decontaminated_counts + out$ambient_counts + out$background_counts, adt_counts, check.attributes = FALSE).
If a decontaminated count is larger, the two matrices are probably not in the same cell or gene order."
DRAFT REPLY: "Thanks both. We checked the implementation: a decontaminated count can't exceed the observed count. Each count is split by an
estimated probability between 0 and 1, and the parts add back to the original exactly. @ivan-pourmir, the explanation in your image doesn't
describe decontX's code: no prior term is added to the native counts. Increases after normalization are expected for each cluster's own
markers, because the cell's total shrinks. The spikier "Decontaminated" curves in plotDensity() are partly a plotting artifact: each curve
gets its own smoothing bandwidth. We'll fix that and add an FAQ entry. If you find a raw count where decontaminated > observed with matching
cell and gene order, please reopen with a small example."

### #43 -> reply now; fix A1-02 + A8-01 in I-b; FAQ III-a
Repro s04: extra genes -> exact error; raw[rownames(x), ] runs; reversed genes silently wrong (max diff 0.59). z = rep(1, n) refused even with
background (A8-01); rep(1:2, length.out = n) + background MAE 0.0067 vs 0.0065 with guard bypassed; estimated z (6 spurious clusters) MAE 0.029.
1.4 GiB sparse->dense warning likely the k-means fallback (A8-02).
DRAFT REPLY: "Thanks, Ethan, and sorry for the cryptic error. background must contain the same genes, in the same order, as the filtered matrix.
decontX matches barcodes (columns) between the two but not genes (rows). The unfiltered matrix usually has more genes than a filtered Seurat
object, so the contamination profile has the wrong number of rows. That's why equalizing the counts didn't help. For now, subset the raw
matrix to your genes first:
  counts <- LayerData(seuratObj, assay = "RNA", layer = "counts")
  raw <- Read10X("raw_feature_bc_matrix")
  raw <- raw[rownames(counts), ]
  result <- decontX(counts, background = raw)
Always subset by rownames(): if the genes are only in a different order, the current version runs without error but gives wrong estimates.
For a single cell type: without background, decontX learns contamination by comparing clusters, so it invents subclusters. That's the odd
result you saw. The 1.4 GiB dense-coercion warning probably comes from its k-means fallback when it finds only one cluster. With background,
the contamination profile comes from the empty droplets, so clusters matter much less. decontX currently refuses a single cluster, but two
arbitrary groups give essentially the same result in our tests: z <- rep(1:2, length.out = ncol(counts)), then decontX(counts, z = z,
background = raw). We plan to check gene names with a clear error, and to allow a single cluster when background is given."

### #45 -> FAQ + reply
Evidence: batches independent (A1-23); 10k ~7.5 s (89% UMAP); 50k 4 batches 36 vs 26 s; peak ~3x dgC; batched ~7.5 copies (A5-05); 100k x4 4.6 GB RSS.
Caveats: >= 2 clusters per batch; 1-cell batch crash (A1-08); char z per-batch codes until A1-09; no column names until Phase II.
DRAFT FAQ (after A1-09): "How do I run decontX on many samples? Contamination differs between samples, so fit each sample separately with batch.
Don't treat merged samples as one: sce <- decontX(sce, z = sce$cell_type, batch = sce$sample). Use raw counts (not integrated, scaled or SCT
values) and keep them sparse (dgCMatrix); never convert the merged matrix to a dense matrix. Put your annotations in z. decontX splits them by
sample, so every sample must contain at least two of the cell types (and more than one cell). The returned decontX_clusters are named
"<batch>-<label>", for example "s1-T cells". With empty droplets, pass background and bgBatch per sample. Runtime grows roughly linearly with
the number of cells: most of it is the PCA/UMAP step, which currently runs even when z is given (about 7-8 s per 10,000 cells on a laptop).
If memory is tight, loop over samples, run decontX on one sample at a time, and cbind() the decontaminated matrices. Each sample is fitted
independently, so the results are the same as with batch." Reply: FAQ + "for Seurat v5 split layers, take each sample's counts layer".

### #46 -> close with comment (fixed in devel 1.11.1: c3d77f7, ee04635, PR #52; not in RELEASE_3_23)
DRAFT CLOSING: "Fixed in decontX 1.11.1 (Bioconductor devel), so it ships in the next Bioconductor release (3.24). decontPro only uses the mean
of the variational approximation, which doesn't depend on the posterior draws, so it now requests the minimum two draws instead of 1000.
Results are unchanged. On a 50 x 400 test set the temporary file dropped from 444 MB to 2 MB and peak memory from 3.0 to 1.2 GB. At your size
the 1000 draws meant roughly 225 GB of CSV, which is gone. rstan still writes a smaller CSV to tempdir(), about 90 bytes per ADT x droplet
(about 1 GB for 215 x 52,222). It is removed when R exits; set TMPDIR before starting R if your temp space is small. We plan to delete it as
soon as decontPro finishes and to expose seed and iter. To use the fix now, install the devel version (BiocManager::install("decontX",
version = "devel")) or the GitHub devel branch."

### #39 -> close (Bioc 3.21 SummarizedExperiment 1.38.0 ExpData class bug; https://stat.ethz.ch/pipermail/bioc-devel/2025-May/020994.html)
DRAFT CLOSING: "This error isn't from decontX. Bioconductor 3.21 had a short-lived bug: SummarizedExperiment 1.38.0 briefly defined an ExpData
class, and SingleCellExperiment builds from that time remembered it. It was fixed in SummarizedExperiment 1.38.1 and a rebuilt
SingleCellExperiment. Update with BiocManager::install() and check with BiocManager::valid(). If it persists, run
BiocManager::install("SingleCellExperiment", type = "source", force = TRUE). Closing; please reopen if the error remains after updating."

### #55 OPEN (test-decon.R:112; I-a/I-c, A3-07). Baseline correction: test warnings = expect_equivalent (#55) + CreateSeuratObject (test-decontPro.R:118).
### #56 FIXED (52de8a0 via #63; man/decontX-package.Rd:6 \value) -> close. Stale RStan citation = A7-12.
### #57 all 4 OPEN: R dep 4.3.0 (DESCRIPTION:59; decision, one line); chunk labels (decontPro.Rmd 10/12 unlabeled, decontX.Rmd L310; III-a);
   decontXcounts.Rd examples (III-a); function length (I-d). Keep open as tracker.

### Findings
A8-01 P2 I-b: one-cluster guard (decon.R:887-893) refuses K=1 even with background (valid model; #43 use case; K=1+bg cor 0.9985 MAE 0.0065);
  without bg K=1 undefined; error doesn't name the batch -> allow per batch when it has background; name batch otherwise. NEWS.
A8-02 P2 I-e/I-b: k-means fallback (decon.R:936-938) on as.matrix(t(normed)) all genes dense -> memory; run on PCA components; changes fallback results (decision); keep legacyInit.
A8-03 P2 III-a: 2026-09 doublet draft order refuted -> use #44 draft.
A8-04 P2 III-a: draft "z = rep(1, ncol(x))" errors today -> rep(1:2, ...) until A8-01.
A8-05 P3 decision: #46 fix not backported to RELEASE_3_23.
A8-06 P3 audit: D7-03 line refs stale (now stan_helpers.R:103-110).
