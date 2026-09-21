# Drafts from the 2026-09 comprehensive review

Companion to `2026-09-comprehensive-review.md`. These are starting
points for the roadmap waves that own them — review before landing.

## 1. NEWS.md backfill (insert between 1.11.1 and 0.99.5)

Ancestry verified against the bioc/RELEASE_3_18..3_23 branches. Note the
DecontPro feature commits shipped in **1.4.1** (a post-release patch on
RELEASE_3_20, 2025-02-16), not 1.4.0.

```markdown
# decontX 1.10.0

* No significant changes (Bioconductor 3.23 release).

# decontX 1.8.0

* No significant changes (Bioconductor 3.22 release).

# decontX 1.6.0

* No significant changes (Bioconductor 3.21 release). Includes the
DecontPro changes released in 1.4.1 on the Bioconductor 3.20 branch.

# decontX 1.4.1 (2025-02-16)

* `decontPro()` gained an `ambient_counts` argument, so an empty-droplet
count matrix can be supplied to estimate the ambient contamination
profile empirically instead of inferring it from the filtered droplets.
* **Breaking:** the first argument of `decontPro()` was renamed from
`object` to `filtered_counts` for consistency with the new
`ambient_counts` argument. Code that passed the count matrix by name
(`decontPro(object = x, ...)`) must be updated; positional calls are
unaffected.
* The estimated ambient contamination profile is now returned in the
`parameters` element of the `decontPro()` output (`p_est`).
* Updated deprecated array syntax in the Stan model, fixed a C compiler
warning in `matrixSums`, and corrected a `deconPro` typo in the
DecontPro vignette.

# decontX 1.4.0

* No significant changes (Bioconductor 3.20 release).

# decontX 1.2.0

* No significant changes (Bioconductor 3.19 release).

# decontX 1.0.0

* First Bioconductor release (Bioconductor 3.18), providing `decontX()`
for ambient RNA decontamination of single-cell RNA-seq counts and
`decontPro()` for single-cell protein expression (CITE-seq/TotalSeq)
data.
```

## 2. inst/CITATION

Both entries Crossref-verified. The DecontPro paper's final form is NAR
52(1):e4 (online 2023-11-17, in-issue January 2024); the repo currently
says "Yin et al. 2023" — if keeping 2023, also update README's link to
the DOI form either way.

```r
citHeader("To cite the decontX package, please cite the method you used.")

bibentry(
  bibtype  = "Article",
  title    = "Decontamination of ambient RNA in single-cell RNA-seq with DecontX",
  author   = c(person("Shiyi", "Yang"),
               person("Sean E.", "Corbett"),
               person("Yusuke", "Koga"),
               person("Zhe", "Wang"),
               person("W. Evan", "Johnson"),
               person("Masanao", "Yajima"),
               person("Joshua D.", "Campbell")),
  journal  = "Genome Biology",
  year     = 2020,
  volume   = 21,
  number   = 1,
  pages    = "57",
  doi      = "10.1186/s13059-020-1950-6",
  header   = "If you used decontX() for single-cell RNA-seq data, cite:"
)

bibentry(
  bibtype  = "Article",
  title    = paste("Characterization and decontamination of background noise",
                   "in droplet-based single-cell protein expression data",
                   "with DecontPro"),
  author   = c(person("Yuan", "Yin"),
               person("Masanao", "Yajima"),
               person("Joshua D.", "Campbell")),
  journal  = "Nucleic Acids Research",
  year     = 2024,
  volume   = 52,
  number   = 1,
  pages    = "e4",
  doi      = "10.1093/nar/gkad1032",
  header   = "If you used decontPro() for single-cell protein expression data, cite:"
)
```

## 3. FAQ sketches (vignette sections)

For `vignettes/decontX.Rmd` (new "Frequently asked questions" section)
and `vignettes/decontPro.Rmd` (after the `plotDensity` chunk):

- **Many samples / memory (#45):** do not `cbind` samples into one dense
  matrix and hope; pass `batch` — decontX subsets and fits each sample
  independently (statistically correct, one dense batch in memory at a
  time). Supply existing annotations via `z`; returned
  `decontX_clusters` are prefixed `"<batch>-<k>"`. If still too large,
  loop per sample and cbind only the decontaminated results
  (equivalent). NOTE: land after the D1-01 numeric-batch fix.
- **`background` requirements (#43):** `x` and `background` must have
  the same genes in the same row order (`background <-
  background[rownames(x), ]`); columns (barcodes) are reconciled
  automatically, rows are not (until the D1-03 validation lands). For a
  single-population dataset, pass `z = rep(1, ncol(x))`.
- **Doublets (#44):** run decontX first on the filtered counts, then
  doublet detection on decontaminated counts; doublets in the input are
  absorbed as high-contamination cells rather than corrupting cluster
  profiles. **Verify the recommended order against Yang et al. 2020
  methods before publishing** (inferred from model structure).
- **Decontaminated values look higher (#40) — CORRECTED per D7-03:**
  decontX/decontPro split each observed count multiplicatively
  (`counts × rate/(sum of rates)`), so a decontaminated count can
  **never** exceed the observed count, and
  `decontaminated + ambient + background == counts` exactly. Taller
  density curves at high expression in `plotDensity` are a per-group
  kernel-bandwidth artifact (the decontaminated negative population is
  compressed → less smoothing); normalized (CLR) values can genuinely
  increase for retained markers — the desired signal-to-noise effect.
  Self-check: `max(decontaminated_counts - counts) <= 0` on aligned
  matrices. Do NOT use the earlier "model expectation, not a
  subtraction" framing — it is wrong for this implementation.
- **Installation errors like "superclass ExpData not defined" (#39):**
  stale mixed library after an R upgrade; run `BiocManager::valid()` and
  reinstall flagged packages.

## 4. Proposed codecov.yml

```yaml
comment:
  layout: "reach, diff, flags, files"
  behavior: default
  require_changes: true

coverage:
  precision: 2
  round: down
  range: "60...90"
  status:
    project:
      default:
        target: 80%          # ratchet up as the test plan lands
        threshold: 1%        # absorbs stochastic-suite jitter
        if_ci_failed: error
    patch:
      default:
        target: 80%
        threshold: 5%
        if_ci_failed: error

ignore:
  - "R/RcppExports.R"
  - "R/stanmodels.R"
  - "src/RcppExports.cpp"
  - "src/stanExports_shrinkage.cc"
  - "src/stanExports_shrinkage.h"
  - "tests/**"
  - "vignettes/**"
  - "man/**"
  - "dev/**"
  - "inst/**"
```

Plus: `Config/Needs/coverage: covr` in DESCRIPTION (metadata, not a
dependency — maintainer's call whether it needs an ADR), a `make
coverage` target, checkout@v4 + upload-artifact@v4, CODECOV_TOKEN +
fail_ci_if_error in test-coverage.yaml, and
`_R_CHECK_FORCE_SUGGESTS_: false` in R-CMD-check.yaml.

## 5. Test plan (18 files, condensed)

Fast = <0.5s, medium = 0.5–3s, slow = >3s. Full assertion detail is in
the audit's dimension-4 section and the review transcript.

| File | Runtime | Guards |
|---|---|---|
| helper-decontX.R | — | fixture constructors (sim_small/sim_sce/adt_counts/fake_vb_fit), custom expectation |
| setup.R | — | suite options, withr teardown |
| test-decon-core.R | medium | recovery oracle vs simulator truth; log-lik non-decreasing; convergence; decontXcounts <= counts |
| test-decon-init.R | medium | existing scrapper-migration blocks (preserve) + seed=NULL, kmeans fallback, varGenes edge |
| test-decon-batch.R | medium | D1-01 numeric labels, 1-cell batch, batch length, per-batch estimate/cell alignment, multi-batch SCE UMAP |
| test-decon-background.R | medium | D1-03 row mismatch, bgBatch guards, eta structure, no-colnames warning |
| test-decon-classes.R | medium | matrix/dgC/DelayedArray/SCE round-trips agree; decontXcounts getter/setter |
| test-decon-validation.R | fast | snapshot every user-facing stop/warning; maxIter=0, iterLogLik=0, bad seed, NA z, gap-level z, label preservation |
| test-simulateContamination.R | fast | identities, delta scalar branch, determinism |
| test-retrieveFeatureIndex.R | fast | full branch coverage (celda_functions 30%→~90%) |
| test-plot-decontx.R | fast | structure/data assertions per no-vdiffr policy; groupClusters; batch=; log1p; labelBars |
| test-plot-decontx-percent.R | fast | .calculateDecontXBarplotPercent exact values; gap-level z; threshold boundary |
| test-plot-decontPro.R | fast | patchwork structure, y-clamp, file= via withr; shared-bandwidth regression (D7-04) |
| test-decontPro.R | 1 slow (~2s) + fast | tiny-fixture identity + visibility; Stan-free guard errors; mocked dispatch incl. Seurat (D3-03) |
| test-matrixSums.R | fast | 3-branch cross-validation vs rowsum oracle; NA-group block (skip until D2-01 guard lands — segfaults) |
| test-matrixNorm.R | fast | fastNormProp == prop.table; equivalence vs celda::normalizeCounts (gates D3-01 swap) |
| test-decontx-kernels.R | fast | all C++ guard messages; estimate_eta/delta contracts; D2-02 min(z) regression; empty-cluster NaN |
| test-utils.R | fast | .logMessages (withr tempfile), .rdirichlet |
```
