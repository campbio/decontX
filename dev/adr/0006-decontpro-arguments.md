# ADR 0006: Rename decontPro's arguments to camelCase and expose its fit controls

## Status

Proposed (2026-10-27). A stub for maintainer review: the decisions below
were made during the Phase 0 review (2026-10-18 to 2026-10-27); the open
points at the end are not yet decided.

## Context

decontPro's arguments are snake_case (`cell_type`, `ambient_counts`,
`delta_sd`, `background_sd`), while decontX and the rest of the package
use camelCase. The review found more:

- `...` swallows unknown arguments, so `seed`, a `delta_SD` typo or
  `batch` are silently ignored (A1-13). `altExpName` is swallowed the
  same way (A4-12).
- The fit's seed, iteration count and number of posterior draws are
  hard-coded and not exposed (D7-01).
- Stan prints about 45 lines to the console, including "EXPERIMENTAL
  ALGORITHM" and "MAY BE DIVERGING" (A7-03). There is no quiet option.
- rstan writes a CSV to `tempdir()`, about 90 bytes per ADT x cell, which
  stays until R exits (A5-08).
- The generic dispatches on both `filtered_counts` and `cell_type`
  (A4-12). With an alias for `cell_type`, a call that uses the alias
  would dispatch on a missing `cellType`.
- The empty-droplet argument `ambient_counts` doesn't match decontX's
  `background`, and its rows are not matched to the counts (A1-12).

The container methods in ADR-0005 change decontPro's return type in the
same devel cycle, so one Bioc release can carry the whole decontPro API
change.

## Decision

We will make camelCase decontPro's primary argument names and keep the
snake_case names as permanent aliases.

**Aliases.** Each existing snake_case argument maps to its camelCase
name (for example `cell_type` to `cellType`, `delta_sd` to `deltaSd`,
`background_sd` to `backgroundSd`). The aliases are permanent, never
deprecated. Giving both names for one argument is an error. Decided
2026-10-18 (area 4, decision 7).

**`background`.** The empty-droplet argument becomes `background`,
matching decontX, with `ambient_counts` as a permanent alias (decided
2026-10-27). `backgroundCounts` and `ambientCounts` were rejected
because they clash with the output elements `background_counts` and
`ambient_counts` (`R/stan_helpers.R:109-115`). The help text separates
`background` (counts from empty droplets) from the model's background
component, whose prior is set by `backgroundSd`.

**Fit controls.** `iter`, `seed` and `output_samples` are exposed
(decided 2026-10-19). A quiet option hides Stan's console output
(decided 2026-10-18; tests silence Stan on their own until then).

**Unknown arguments.** An unknown argument in `...` is an error, not
silently ignored (A1-13).

**Generic signature.** `setGeneric("decontPro", ..., signature =
"filtered_counts")`, so dispatch uses the first argument only and the
`cell_type` alias works (decided 2026-10-27). This changes no
behaviour.

**Temporary file.** rstan's CSV is removed as soon as the fit ends
(A5-08).

**No model change.** The Stan model and the VB defaults stay as they
are. Runtime and the iteration cap are documented instead.

## Consequences

- decontPro reads like decontX. Existing scripts keep working through
  the aliases.
- Typos and unsupported arguments fail loudly. Scripts that passed
  ignored arguments now error; NEWS says so.
- Users can reproduce a fit with `seed` and bound its runtime with
  `iter`.
- Temporary disk use no longer grows across fits in one session.
- The help pages document two names for each renamed argument.
- The vignette (Phase III-a) uses the camelCase names, a short
  `background` example (`eval = FALSE`), and the quiet option once it
  exists.

## Open points

- Whether the first argument, `filtered_counts`, also gets a camelCase
  name. It is the dispatch argument, so S4 can't treat it as an alias:
  if it is renamed, `decontPro(filtered_counts = x)` reaches the
  generic only through `...`, and needs a `"missing"` method that
  redirects. Keeping `filtered_counts` as the one snake_case name
  avoids that; most calls pass it by position anyway.
- The name of the quiet option (for example `verbose` or `quiet`) and
  its default.
- The camelCase spelling of `iter`, `seed` and `output_samples` (for
  example `outputSamples`), and their defaults, which must keep today's
  results.
- How the aliases are implemented (argument handling in each method or a
  shared helper), including whether a `"missing"` method is needed for
  the dispatch argument.
- Whether `ambient_counts` keeps its current row-sum semantics under the
  new name, and how its rows are matched (A1-12, fixed in Phase I-b).
