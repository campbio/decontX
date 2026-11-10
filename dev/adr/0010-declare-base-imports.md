# ADR 0010: Declare stats, grid and utils in Imports

## Status

Accepted (2026-11-10)

## Context

decontX calls functions from three packages that ship with R but aren't
attached by every session in the same way: `stats` (for example
`stats::kmeans`, `stats::rbeta`, `stats::rmultinom`), `grid`
(`grid::unit` in the marker percentage plot) and `utils`
(`utils::type.convert` in `.meltMatrix()`). The calls already use `::`,
but DESCRIPTION doesn't list the three packages, so the dependency is
only implied. R CMD check doesn't flag it, so nothing caught it.

The 2026-10 deep review (A4-16 in `dev/audits/2026-10-deep-review.md`)
recorded it, and the maintainer agreed to declare them (open question 2
of the audit; issue #116).

## Decision

We will list `grid`, `stats` and `utils` in Imports, without version
requirements, as they come with R. The calls keep their `::` prefixes,
so NAMESPACE doesn't change and nothing is imported into decontX's
namespace.

Alternative considered: leaving them undeclared. Base R always installs
them, so nothing breaks today, but every other package decontX calls is
declared, and tools that read DESCRIPTION (dependency graphs, renv,
BiocCheck) see an incomplete list.

## Consequences

No change for users: the three packages come with every R installation,
and results, plots and the exported API are unchanged. DESCRIPTION now
lists every package decontX's code calls. celda and singleCellTK are
unaffected.
