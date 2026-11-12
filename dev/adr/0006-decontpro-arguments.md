# ADR 0006: Rename decontPro's arguments to camelCase and expose its fit controls

## Status

Accepted (2026-11-11). Drafted as a stub on 2026-10-27 from decisions
made during the Phase 0 review (2026-10-18 to 2026-10-27); the open
points were settled on 2026-11-11 (see "Resolved points").

## Context

decontPro's arguments are snake_case (`cell_type`, `ambient_counts`,
`delta_sd`, `background_sd`), while decontX and the rest of the package
use camelCase. The review found more:

- `...` swallows unknown arguments, so `seed`, a `delta_SD` typo or
  `batch` are silently ignored (A1-13). `altExpName` is swallowed the
  same way (A4-12).
- The fit's seed, iteration count and number of posterior draws are
  hard-coded and not exposed (D7-01). The draw count was the workaround
  asked for in #46; #52 then set it to 2 and kept draws for `tau_a` only.
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

**Fit controls.** `iter` and `seed` are exposed (decided 2026-10-19).
A quiet option hides Stan's console output (decided 2026-10-18).
`output_samples` is not exposed (decided 2026-11-11, reversing the
2026-10-19 decision): decontPro reads only the mean of the variational
approximation, so the number of draws never changes its results, and it
keeps no draws to return. More draws would only bring back the memory and
disk use of #46, and fewer than 2 crash rstan. Returning posterior draws
would be a separate feature with its own design.

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
- Users can choose the seed of a fit with `seed` and bound its runtime
  with `iter`.
- Temporary disk use no longer grows across fits in one session.
- The help pages document two names for each renamed argument.
- The vignette (Phase III-a) uses the camelCase names, a short
  `background` example (`eval = FALSE`), and the quiet option once it
  exists.

## Resolved points

Settled with the maintainer on 2026-11-11.

- **First argument.** `filtered_counts` keeps its name, the one
  snake_case argument. It is the dispatch argument, so an alias would
  need a `"missing"` method; most calls pass it by position.
- **Quiet option.** `verbose = TRUE`, as in decontX, so the default
  prints what it prints today. With `verbose = FALSE`, Stan's output is
  captured, and shown with `message()` before the error if the fit fails
  (an error message is cut at 1000 characters).
- **Spelling and defaults.** `iter = 50000` and `seed = 12345`, today's
  hard-coded values, so default results are unchanged. They and
  `verbose` come after `...`, so they must be named and a stray
  positional value still errors.
- **Aliases.** The generic is `function(filtered_counts, cellType, ...)`
  with `signature = "filtered_counts"`. The methods keep today's
  positional order (`cellType`, `deltaSd`, `backgroundSd`, `background`),
  and one shared helper reads the aliases from `...`. No `"missing"`
  method is needed.
- **`background` semantics.** Unchanged: its row sums, normalized, give
  the ambient profile, and its rows are matched to the counts as since
  #89 (A1-12).
- **Out of scope.** The output list keeps its names
  (`parameters$delta_sd` and the others); ADR-0005 (#82) reworks the
  return value.
