# Phase I-a regression net (#66)

Status: proposed, 2026-10-30. Roadmap Phase I-a; audit findings A3-01,
A3-05, A3-06, A3-10, A3-16 (in part), A3-17, A3-19, and the "P" rows of
the I-a pinning plan in `dev/audits/2026-10-review-drafts.md`.

## Why

Today the tests only compare results within one run, so nothing shows
when a fix or refactor changes results (A3-01). Before any Phase I-b fix,
this PR stores reference outputs and pins today's behaviour wherever a
later fix will touch it. Every new test passes on today's code; the
regression tests that fail now belong to each I-b fix, written test-first
in that fix's PR.

## What changes

### 1. Committed reference outputs (A3-01, A3-17)

- `tests/testthat/fixtures/make-refs.R` regenerates
  `tests/testthat/fixtures/decontx-refs.rds` (xz). It runs on
  `simulateContamination(C = 100, G = 60, K = 3, seed = 1)` with the
  simulated `z` and `verbose = FALSE`, over five paths:
  - `single`: matrix, unbatched
  - `batched`: `batch` = two halves, "A" and "B"
  - `background`: `background` = a second simulation (`C = 40`,
    `seed = 2`), with the same gene names and different cell names
  - `bgBatched`: batched, plus `background` with a matching `bgBatch`
  - `fixedDelta`: `estimateDelta = FALSE`
- **What each path stores:** `contamination`, `z`, `decontXcounts`
  (dgCMatrix), and, per batch, `theta`, `phi`, `eta`, `delta`,
  `logLikelihood` and `iteration`.
- **What it leaves out:** the UMAP. With `z` supplied, the UMAP doesn't
  feed the EM, and it depends on the platform.
- **Also stored:**
  - the full output of `simulateContamination()` (the default call and
    the `seed = 1` call), to pin the simulator before the I-d
    `rowsum()` swap
  - `.calculateDecontXBarplotPercent()` on the `single` path with labels
    1..K (the A1-07 pin)
- `tests/testthat/test-references.R` reruns each path and compares it
  with `expect_equal(tolerance = 1e-10)`. The 1e-10 allows for libm
  differences across the Bioc builders. Integer parts (`z`, `iteration`,
  simulated counts) must match exactly. All five paths are bit-identical
  on rerun locally, and the EM has no thread-dependent reductions.
- **About 70 KB** without the UMAPs.
- **When a later PR changes results on purpose**, it reruns
  `make-refs.R` and says in the PR which references changed and why.
  The test failure itself shows which element moved.

### 2. Cross-commit identity script (A3-01)

- `dev/compare-commits.R <old-ref> [<new-ref>]` (new ref default
  `HEAD`). It runs locally only, never in CI.
- **How it works:** it installs each commit from `git archive` into its
  own temporary library, runs a fixed set of paths in a fresh R process
  per library (callr), and compares them with `identical()`. It prints
  one line per path: identical, or the largest absolute difference.
- **The paths:**
  - the five above
  - no `z`, for both the matrix and SCE (with the UMAP)
  - `legacyInit = TRUE` (when scater is installed)
  - both `simulateContamination()` calls
  - decontPro on the 10 x 8 fixture, with and without `ambient_counts`
    (real Stan fits)
- It carries the reproducibility checks that can't be committed at
  1e-10: no `z`, legacyInit, and decontPro (maintainer decisions 2 and 9).
- **Cost:** installing a commit compiles Stan, so expect several minutes
  per commit.

### 3. Pins on today's behaviour (all pass now)

- **A3-05:** the numeric batch-label test uses two simulations with
  different seeds, so a batch misalignment would show.
- **A3-06:** each bare expectation gets the pattern of today's message:
  - `test-decon.R:121`: the `simulateContamination()` warning
    ("Only 10 clusters are simulated")
  - `test-decon.R:122`: the `numMarkers` x `K` > `G` error
  - `test-decon.R:284` and `:288`: the `delta` error
  - `test-matrixSums.R:13`, `:14`, `:24` and `:25`: the factor and
    length errors
- **A1-03:** test the two existing `bgBatch` stops: `bgBatch` given
  without `batch` (`decon.R:424`), and `batch` plus `background` without
  `bgBatch` (`:437`).
- **A1-10/A1-17:** test the `varGenes` and `dbscanEps` stops
  (`decon.R:1051`, `:1063`). Both fire even when `z` is given.
- **A3-10 / A1-12:** cover the decontPro `ambient_counts` path for the
  matrix and SCE methods, with `.call_stan_vb` mocked to capture the
  Stan data. This pins `dat$p` as `rowSums(ambient) / sum(ambient)`,
  matched by position as today. No Stan run.
- **A3-19 / A1-06 / A1-14:**
  - dgCMatrix input gives results identical to dense input.
  - DelayedArray(dgCMatrix) input gives identical results and returns a
    DelayedMatrix.
  - Dense input returns a dgCMatrix, which is today's behaviour; #71
    will change this.
- **A3-16** (only the part needed here): new tests keep their state
  local with withr (`local_seed`, `local_tempfile`). Message wording for
  the two misspelled messages that #74 will fix ("bacth", "numebr") goes
  in `expect_snapshot()`; the regexp stays the primary check (decision
  8). On Bioc, snapshots are skipped.

## Not in this PR

- **The HDF5 input test (A3-19):** HDF5Array isn't a declared
  dependency, and R CMD check flags undeclared packages used in tests.
  It goes in when HDF5Array joins Suggests (maintainer decision 6,
  ADR-0005).
- **Regression tests that fail today:** these go in their I-b fix PRs
  (#67 to #74).
- **Moving top-level fixtures, fixing the mutated fixture, the hollow
  tests, renaming the legacyInit test, and the test warnings:** all in
  I-c (#75).
- **plotDensity pins (A1-24):** these go with the #73 plot fixes.

## Docs, NEWS, version

- No user-facing change, so no NEWS entry.
- Version: z + 1. It becomes 1.11.4 if #86 (1.11.3) merges first;
  otherwise 1.11.3. Set when rebasing on devel before the PR.
- No package code, roxygen or vignette changes.

## How it is checked

- Each new test passes on today's code. To show the reference test can
  fail, I'll perturb one path by hand once during development (for
  example `delta = c(10, 11)`), confirm it fails, and not commit that.
- `make test` stays fast: the new tests add about 1 s.
- Coverage doesn't drop.
- `make check-full` and `make bioccheck` before the PR.
- `dev/compare-commits.R devel HEAD` reports every path as identical,
  since no package code changed.
