# ADR 0003: Remove the celda dependency and invert the celda/decontX relationship

## Status

Proposed

## Context

decontX imports celda (DESCRIPTION Imports) but uses exactly **one**
celda function, in exactly **one** place:

```r
# R/decon.R, inside .simulateContaminatedMatrix()
eta <- celda::normalizeCounts(counts = nGByK, normalize = "proportion")
```

`normalize = "proportion"` is a column-wise proportion normalization
(each column divided by its sum) — the same operation the package's own
kept `fastNormProp` C++ routine already performs. Pulling in all of celda
for this is a large, brittle dependency for a one-line need.

Two forces make this more than a routine dependency-diet item:

1. **Deprecation churn.** celda's `normalizeCounts` sits on top of
   scuttle/scater helpers that Bioconductor has been deprecating
   (this already surfaced as `R CMD check` warnings in decontX's
   examples). Each celda-internal change can ripple into decontX even
   though decontX's own use is trivial. (Note: we are *not* dropping
   celda as a reaction to that deprecation — that fix belongs upstream
   in celda and is being handled there. This ADR is a deliberate,
   planned removal, tracked separately from the deprecation.)

2. **A masked duplicate and a dependency-direction problem.** celda
   currently ships a ~5-year-old copy of decontX's functions under the
   same names, which masks the real decontX exports for users who load
   both (see the celda/decontX masking issue). The right long-term
   shape is the **inverse** of today's dependency: celda's DecontX /
   DecontPro entry points should call into (or re-export from) this
   package, and this package should not depend on celda at all. celda
   is being updated in parallel to make that possible.

The hard constraint is **Bioconductor safety**: at no single point may
either package fail for users, and no release may introduce a masking
conflict, a circular Imports edge, or a version-skew breakage. Both
packages release on the same Bioconductor cadence, so the change has to
be staged so that every intermediate state (each package at each release)
is self-consistent.

## Decision

We will remove decontX's dependency on celda, in two coordinated stages.

**Stage 1 — drop celda from decontX (this package, self-contained).**
Replace the single `celda::normalizeCounts(..., normalize =
"proportion")` call in `.simulateContaminatedMatrix()` with a small
local proportion-normalization helper (base R, or a thin wrapper over the
existing `fastNormProp`), verified to produce identical output on the
simulation path, and remove `celda` from DESCRIPTION Imports. This stage
depends on nothing in celda and can land as soon as the replacement is
test-covered. It does **not** by itself resolve the masking issue — celda
still ships its duplicate — but it removes the decontX → celda edge and
the deprecation-ripple exposure.

**Stage 2 — invert the relationship (cross-package, coordinated).**
Work with the parallel celda update so celda's DecontX/DecontPro surface
references this package instead of carrying a stale duplicate, retiring
the name masking. This stage is sequenced across both packages'
Bioconductor release branches so that no release of either package, at
any point, (a) fails for users, (b) masks the other's exports, or
(c) forms a circular Imports dependency. The concrete ordering (which
package changes in which Bioc cycle) is agreed on the celda side and
recorded here when settled.

Alternatives considered:

- **Keep celda for the one call.** Rejected: a whole imported package for
  a one-line column normalization, plus continued exposure to
  celda-internal deprecation churn and the masking problem.
- **Vendor a copy of celda's `normalizeCounts`.** Rejected: it would be
  more code than the operation needs; the "proportion" mode is a
  one-liner we already implement in `fastNormProp`.
- **Do Stage 2 first (invert before dropping the Imports edge).**
  Rejected as the opening move: Stage 1 is safe and independent, so it
  should not wait on cross-package coordination.

## Consequences

- **Easier:** decontX's dependency graph shrinks (celda and its
  transitive dependencies leave Imports); check/build times drop; a
  class of celda-internal deprecation warnings can no longer reach
  decontX's checks. `R/celda_functions.R` (helpers mirrored from celda)
  stays as-is for now and is unaffected by Stage 1.
- **Harder / to watch:** Stage 2 is cross-package and cannot be completed
  from this repo alone — it is gated on the parallel celda work and must
  be sequenced against Bioconductor releases so no intermediate state
  breaks. Until Stage 2 lands, the celda-side masking of decontX exports
  persists for users who load both packages.
- **Follow-up:** implement Stage 1 (local helper + drop Imports + test)
  under this ADR; track Stage 2 as a cross-repo coordination item and
  amend this ADR with the agreed release ordering once fixed. Supersede
  or mark this ADR Accepted when Stage 1 lands.
