# ADR 0002: Migrate cluster initialization from scater/scuttle to scrapper

## Status

Accepted (2026-09-21)

## Context

The Bioconductor 3.24 devel build report for decontX shows four
`'normalizeCounts' is deprecated` warnings from the examples. scuttle
deprecated `normalizeCounts()` and `logNormCounts()` in Bioc 3.24 as
part of the ecosystem-wide migration to the scrapper package (libscran
C++ implementations); deprecated functions may become defunct at the
next release, which would turn these warnings into build errors.
decontX called `scater::logNormCounts()` (a scuttle re-export) in
`.decontxInitializeZ()`, the pipeline that log-normalizes, selects
variable genes, and runs PCA/UMAP to derive initial cell clusters when
`decontX()` is called without `z`. scrapper is available on Bioconductor
release (3.23) with a lighter dependency footprint than scater.

## Decision

We migrate the whole initialization pipeline (normalization → variable
gene selection → PCA → UMAP) to scrapper and make it the default:
`centerSizeFactors` + `normalizeCounts`, `modelGeneVariances` +
`chooseHighlyVariableGenes`, `runPca` (50 PCs, matching scater's
internal default), and `runUmap` (explicitly seeded — deterministic
without R's RNG). scrapper is added to Imports; scater moves from
Imports to Suggests (still used by vignettes and the legacy path).

Because the new pipeline produces a slightly different embedding (HVG
selection by variance-model residuals instead of raw variance; umappp
instead of uwot), a `legacyInit = TRUE` argument keeps the original
scater/scuttle path available so published analyses can be reproduced.
The legacy path loads scater on demand and may emit the upstream
deprecation warnings by design. It is a permanent backwards-compatibility
option, maintained for as long as scater/scuttle continue to provide the
underlying functions (its lifetime is therefore bounded by upstream, not
by us).

Alternatives considered: (a) inlining the log-normalization formula and
keeping `scater::calculateUMAP` — no new dependency and bit-identical
results, but leaves decontX on the deprecated side of the ecosystem
migration and keeps heavy scater in Imports; (b) suppressing the
warnings — rejected, breaks when the functions go defunct.

## Consequences

- `decontX()` results without `z` change slightly across versions
  (documented in NEWS; before/after concordance checked on simulated
  data in the migration PR). Users passing `z` are unaffected.
- Install footprint shrinks: scater and its dependency tree leave the
  hard requirements; scrapper brings beachmat/BiocNeighbors/SparseArray.
- The UMAP step becomes deterministic on the default path even when
  `seed = NULL`.
- `legacyInit` is kept indefinitely for backwards compatibility. If
  scuttle ever makes the deprecated functions defunct, the legacy path
  will fail with upstream's defunct error — an accepted risk outside our
  control; the option is only removed if that happens.
- The `runUmap` seed-argument API differs between scrapper <= 1.4
  (`seed`) and >= 1.5 (`initialize.seed`/`optimize.seed`); the code
  supports both via a formals check until old versions age out.
