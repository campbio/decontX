# Validate decontX inputs (#67)

Status: approved by the maintainer, 2026-10-30. Roadmap Phase I-b; audit
findings A1-01, A1-10, A1-16, A1-17, A1-18, A1-19, A5-07, A6-01 and A7-13.

## Why

Several inputs aren't checked. Some give silently wrong results: a `z`
shorter than the number of cells is padded with NA, which becomes its own
cluster. Others give cryptic errors from deep inside R, withr or
scrapper. All of these were reproduced on devel at cf862bb; the
correlation with the true contamination was 0.999 for a valid `z`.

- **Wrong results, silently:**
  - `z` of length 90: correlation 0.51
  - `z = "cluster"`: correlation −0.02
  - `z` with NAs: correlation 0.64
  - `convergence = "x"`: stops after 1 iteration
  - `delta = c(0, 0)`: correlation 0.009
  - negative counts: NaN in the output
- **Invalid values accepted silently:** `maxIter = "a"`,
  `dbscanEps = "a"`, and logical flags set to strings or NA.
- **Cryptic or misleading errors:**
  - `batch` of the wrong length, or with NA
  - `maxIter = 0`, `iterLogLik = 0`
  - NA in `varGenes`, `dbscanEps` or `delta`
  - `seed = "a"`
  - an `assayName` that isn't an assay of `x`

**Root cause:** validation is scattered and late.
- `z` is checked only per batch, after `z[batch == bat]` has already
  padded it.
- `batch` and the scalar arguments are never checked.
- `.processvarGenes`, `.processdbscanEps` and `.checkDelta` let NA, or 0
  for `delta`, through.

## What changes

1. **`.checkDecontXArgs()`.** A new internal function that both methods
   call first, before background processing or any other work. Each
   message names the argument and the problem. It checks:
   - `z` and `batch`: NULL, or one value per cell with no NA
   - `maxIter`: a single whole number, at least 1, or `Inf` (run until
     convergence, which worked before)
   - `iterLogLik`: a single whole number, at least 1, or `Inf` (log
     likelihood only at convergence, which worked before)
   - `convergence`: a single non-negative number, above 0 when
     `maxIter` is `Inf` (otherwise the EM never stops)
   - `estimateDelta`, `legacyInit` and `verbose`: `TRUE` or `FALSE`
   - `seed`: NULL or a single number in R's integer range
   - `varGenes`: a single whole number, at least 2, or `Inf` (all genes,
     which worked before)
   - `dbscanEps`: a single non-negative number
   - `delta`: two finite, positive numbers
   - `logfile`: NULL or a single character string

   The `varGenes`, `dbscanEps` and `delta` messages keep their current
   wording, so the existing pins keep matching.
2. **Counts:** `.checkCountsDecon()` errors on negative counts and uses
   `anyNA()` (Matrix's method) instead of `sum(is.na())` (A5-07). It runs
   once on the whole counts matrix, before any batch, and on the
   background counts.
3. **Assay names:** the SCE method checks `assayName` against
   `assayNames(x)`, and `bgAssayName` against the background's assays,
   before anything runs (A6-01). An assay index still works. When
   `bgAssayName` defaults to `assayName`, the message names
   `assayName`.
4. **Docs (A7-13):**
   - `varGenes` and `legacyInit` also build the UMAP stored in the
     results when `z` is given.
   - `dbscanEps` is used only when `z` is NULL.
   - The zero-count error no longer tells users to supply `z`, since
     initialization runs either way until #79 (A1-16).

## Tests

- One test per check, written first, that fails on its assertion on
  devel.
- No change for valid input: the stored references still match, and
  `dev/compare-commits.R devel HEAD` reports every path identical.

## NEWS, version, release

- NEWS entry: invalid inputs now stop with an error that names the
  argument.
- Version 1.11.5.
- Devel only. The wrong results come only from malformed input.

## Effect on related packages

celda forwards its arguments to `decontX()`. Calls with valid values are
unaffected; invalid values that used to run now stop with an error.

## Out of scope

- background row matching and `bgBatch` (#68)
- the class of the returned `z` (#70)
- dgR, dgT and data.frame input (#71)
- consolidating parameters into one constructor (ADR-0004, I-d)
