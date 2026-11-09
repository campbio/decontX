# ADR 0008: Drop plyr and reshape2 from Imports

## Status

Accepted (2026-11-08)

## Context

decontX imported plyr for one function, `mapvalues()`, used twice: to
number the cluster labels in `.processCellLabels()` and to renumber the
simulated clusters in `simulateContamination()`. It imported reshape2 for
`melt()`, used at five sites in the plotting code. Both packages are
retired: their authors recommend dplyr and tidyr, and they only get
fixes needed to stay on CRAN. Each one is an install dependency for every
decontX user, for a few lines of code that base R can do.

The 2026-10 deep review (A4-09, A4-10 in
`dev/audits/2026-10-deep-review.md`) found that the obvious replacement
for `mapvalues()`, `match(z, unique(z))`, is not identical.
`mapvalues()` keeps the type of the labels, so character labels get
character codes ("1", "2", ..., "10"), which `as.factor()` sorts as text
("1", "10", "2"). The EM numbers the clusters by these levels, so with
10 or more character labels a plain `match()` would reorder the columns
of phi and eta and change contamination in the last digits.

## Decision

We will remove plyr and reshape2 from Imports and use base R:

- Labels: `z[] <- match(z, unique(z))`, which keeps the type of `z`
  exactly as `mapvalues()` did, so the cluster order and every result
  stay bit-identical. In `simulateContamination()`, where `z` is always
  an integer, `match(z, unique(z))` alone.
- Matrices: an internal `.meltMatrix()` that reproduces
  `reshape2::melt()` on a matrix, including its conversion of dimnames
  with `utils::type.convert()` (numbers stay numbers, other names become
  factors in their order, missing names become indices). The marker
  plots depend on these column types.
- The decontPro plots: an internal `.meltConDecon()` that stacks the
  counts before and after decontamination, as `melt()` did on their
  two-column data frame.

Alternatives considered: tidyr (`pivot_longer()`) and data.table
(`melt()`) would replace one dependency with another, and neither
converts dimnames the same way. A simpler matrix melt with
`as.data.frame(as.table(m))` turns every dimname into a factor, which
changes the plot data for numeric-like cell names.

## Consequences

Two fewer Imports, with no change in results, numbers or plots: the
regression references are identical. The two helpers are about 30 lines
of code we now maintain, with tests that pin their output. Nothing
exported changes, so celda and singleCellTK are unaffected. Users who
loaded plyr or reshape2 only through decontX now have to install them
themselves.
