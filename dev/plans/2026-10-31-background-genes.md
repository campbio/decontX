# Match background rows by name, and check `bgBatch` (#68)

Status: approved by the maintainer, 2026-10-31. Roadmap Phase I-b; audit
findings A1-02, A1-03 and A1-12. Refs #43.

## Why

Reproduced on devel at 0040d2c. decontX: 200 simulated cells, 300 genes,
100 empty droplets; correlation with the true contamination is 0.994
when the background rows match.

- **decontX `background` (A1-02):**
  - reversed rows run: correlation 0.875, contamination off by up to 0.45
  - shuffled rows run: correlation 0.847, off by up to 0.52
  - one extra or one missing gene: the C++ error "The number of rows in
    'eta' must be equal to the number of rows in 'counts'" (#43)
  - row names that don't overlap at all run
- **`bgBatch` (A1-03), batches A and B (correlation 0.992 when matched):**
  - `bgBatch` of C and D runs and ignores the background: 0.856
  - `bgBatch` of only A leaves batch B without background: 0.663
  - a `bgBatch` of length 2 is recycled without a warning
  - NA gives "missing value where TRUE/FALSE needed"
- **decontPro `ambient_counts` (A1-12), 10 x 8 fixture, real fit:**
  reversed rows run and change the decontaminated counts by up to 0.53;
  ADT names that don't overlap run.

**Root cause:** the contamination profile is `rowSums()` of the
background by position (`decon.R:809`, `decontPro.R:137-147`); row names
are never read. Background columns are split with
`batchBackground == bat` (`decon.R:522`), which gives a 0-column
background for a batch with no match (a flat profile) and recycles short
vectors. The `bgBatch` length check in `.checkBackground()` runs only
when duplicate barcodes were removed.

## What changes

1. **`.matchBackgroundRows()`**, shared by decontX (`background`) and
   decontPro (`ambient_counts`):
   - Both have row names:
     - identical names: used as they are (this also covers duplicated
       gene names);
     - otherwise, duplicated names in either: error, since matching by
       name would be ambiguous;
     - names of the counts missing from the background: error naming the
       number missing and the first few;
     - extra background rows: dropped, with a log message giving the
       number dropped and the first few names;
     - rows reordered to match the counts.
   - Either has no row names: matched by position when the row counts
     are equal (today's behaviour); otherwise an error. decontPro keeps
     "different number of ADTs" in its message.
   - decontX calls it in `.decontX()` before any work; decontPro in
     `.decontPro()`, in place of the row-count check. decontPro logs with
     `message()`.
2. **`bgBatch`:**
   - `.checkBackground()` first checks that `bgBatch`, when given, has
     one value per background column and no NA.
   - `.decontX()` checks that every `bgBatch` value is a `batch` value,
     and that every batch has background droplets. This runs after
     barcodes shared with `x` are removed, so it also catches a batch
     that loses all its droplets that way.
3. **Empty background:** an error when no droplets are left after
   removing barcodes found in `x`. Today the profile silently becomes
   flat.
4. **Docs:** the `background`, `bgBatch` and `ambient_counts` params say
   rows are matched by name, and that every batch needs background.

## Tests

Written first; each fails on its assertion on devel.

- decontX: reversed rows, and extra genes, give results identical to
  matched rows; matrix and SCE methods.
- One test per new error: missing genes, no overlap, no names with a
  different row count, ambiguous duplicated names, `bgBatch` length, NA,
  a value not in `batch`, a batch without background, empty background.
- The extra-gene log message.
- decontPro, with Stan mocked: reversed `ambient_counts` rows give the
  same profile as matched rows; missing ADTs error. The existing
  "different number of ADTs" test keeps passing.
- No change for valid input: the stored references (same row names, same
  order) still match, and `dev/compare-commits.R devel HEAD` reports
  every path identical.

## NEWS, version, release

- NEWS entry: background rows are matched by gene name, `bgBatch` is
  checked, and reordered rows now give correct results.
- Version 1.11.6.
- Devel only, with the other Phase I-b fixes.

## Effect on related packages

celda forwards `background` and `bgBatch` to `decontX()`. Inputs whose
rows already match give identical results; reordered rows now give
correct results; non-matching names or `bgBatch` values now stop with an
error.

## Out of scope

- one cluster with background, and the k-means fallback (#69)
- a matrix `background` passed to the SCE method
- decontPro argument names (`background` alias, camelCase; #81)
