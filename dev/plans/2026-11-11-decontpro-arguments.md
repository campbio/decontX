# decontPro arguments: camelCase, `background`, fit controls (#81)

Approved in chat 2026-11-11, with ADR-0006 (now Accepted). The maintainer
dropped `output_samples` from the API.

## Change

- **Names.** The methods take `filtered_counts, cellType, deltaSd = 2e-5,
  backgroundSd = 2e-6, background = NULL, ...`, so positional calls keep
  working. `filtered_counts` keeps its name. `cell_type`, `delta_sd`,
  `background_sd` and `ambient_counts` are permanent aliases, read from
  `...` by one shared helper. Giving an argument and its alias (or an
  alias twice) is an error; any other name in `...` is still an error.
- **Fit controls.** `iter = 50000`, `seed = 12345` and `verbose = TRUE`
  come after `...`, so they must be named and a stray positional value
  still errors. The defaults are today's hard-coded values, so results
  don't change. `iter` is a whole number of at least 1, `seed` a whole
  number from 0 to `.Machine$integer.max`, `verbose` TRUE or FALSE.
  `iter` is at most `.Machine$integer.max` too. Unknown arguments are
  reported before anything in `...` is evaluated.
- **Quiet option.** `verbose = FALSE` captures Stan's console output.
  When the fit fails, the captured output is shown with `message()`
  before the error (an error message is cut at 1000 characters).
- **Generic.** `setGeneric("decontPro", function(filtered_counts,
  cellType, ...), signature = "filtered_counts")`.
- **Temporary file.** `.callStanVb()` passes its own `sample_file` to
  `rstan::vb()` and deletes it when the fit ends.
- **Docs.** Each renamed argument documents its alias. `deltaSd` and
  `backgroundSd` are documented as the SDs of the half-normal priors on
  the spread of the ambient and background levels, not variances.
  `background` is described as empty-droplet counts, separate from the
  model's background component. Error messages use the camelCase names.
- **Unchanged.** The Stan model, the VB settings, the output list and its
  names (`parameters$delta_sd` and so on; #82 reworks the return value),
  the number of draws (2), and `background`'s row-sum semantics. The
  decontPro vignette still runs through the aliases; #83 part 2 rewrites
  it.

## Tests (written first)

- Each alias gives the same call into the core as its camelCase name
  (core mocked), for the matrix, SCE and Seurat methods.
- Both names, or an alias twice, is an error naming both.
- Unknown arguments still error; existing tests keep passing.
- `iter` and `seed` reach `rstan::vb()` (mocked); invalid `iter`, `seed`
  and `verbose` error before any fit.
- `verbose = FALSE` prints nothing on a real fit (10 x 8 fixture); a
  failed quiet fit shows Stan's output in the error.
- rstan's CSV no longer exists after a fit.

## Version, NEWS, checks

1.11.33, with a NEWS entry. `dev/compare-commits.R` should show
identical results on all paths. `make test`, `make coverage`,
`make check-full`, `make bioccheck`.
