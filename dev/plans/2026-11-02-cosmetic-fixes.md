# Cosmetic fixes (#74)

Approved in chat, 2026-11-02. Audit findings A1-21 and A1-26
(`dev/audits/2026-10-deep-review.md`). Version 1.11.12.

## Changes

1. `simulateContamination()`, `NRange`: draw each cell's total from the
   range with `sample.int()` on its length. When min < max the draws are
   the same as before, so the default simulation and the stored references
   don't change. When min equals max, every cell gets that total (it used
   to draw from `1:min` through `sample()`'s length-1 behavior).
2. `simulateContamination()`, `K`: stop early when `K < 2`, since the
   simulated contamination comes from the other populations (it crashed
   with "Division by 0"). Stop with a clear message when sampling leaves
   only one population (for example `C = 1`).
3. `retrieveFeatureIndex()`: a feature that is `NA` never matches (it went
   to a row named `NA`); it goes into the "not present" warning. `by` must
   be a single non-NA string, checked first.
4. `decontX()` `runParams`: add `iterLogLik` (after `convergence`) and
   `seed` (after `dbscanEps`), per the 2026-10-18 decision.
5. Typos: messages ("bacth", "backgroud", "numebr") and roxygen
   ("lenght", "paramters", "percentation of contamation", "retived",
   "retreived", "retrived"). The `simulateContamination()` `@return` names
   `nativeCounts`/`observedCounts` (it said `nativeMatirx`/
   `observedMatrix`). The `bgBatch`-without-`batch` error is reworded to
   name both arguments. `make docs` regenerates the help pages.

## Tests

`test-decon.R`, each failing first: equal `NRange` values; `K = 1` and
`C = 1` errors; `NA` features with both `exactMatch` settings; bad `by`;
`runParams` holding `seed` and `iterLogLik`. The snapshots that recorded
the misspelled messages are updated, and the regexp on the reworded
`bgBatch` error follows the new wording.

## Results

decontX results are unchanged (check with `dev/compare-commits.R`);
`runParams` gains two fields, which celda's forwarding layer passes
through. Simulations change only for one-value `NRange` and `K = 1`.

## Related packages

celda ships an identical copy of `retrieveFeatureIndex()`; the `NA` and
`by` fixes need a matching celda issue (drafted for review, not filed).
