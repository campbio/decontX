# `.process_stan_vb_out()` tidy-up (#80, part 3)

Approved in chat 2026-11-08. Audit finding A4-20.

## Problem

`.process_stan_vb_out()` in `R/stan_helpers.R` splits decontPro's counts
into native, ambient and background parts, three N x M matrices (ADTs x
droplets). It:

1. computes the denominator `ambient_rate_est + cell_rate_est +
   background_rate_est` three times, two extra full-size matrices and four
   extra element-wise additions;
2. assigns `background_est` from `val$background` twice; the second
   assignment changes nothing;
3. reuses `unscaled_rates` for the ambient profile and then the native
   profile.

## Change

- Compute the denominator once (`total_rate_est`) and divide by it in all
  three parts. Each line keeps the form `counts * rate / total_rate_est`,
  so R evaluates it in the same order as before.
- Drop the second `background_est` assignment.
- Name the two profiles `ambient_profile` and `native_profile`.
- No change to the formula, the outputs or their names.
- Out of scope: the 0/0 (NaN) case when all three rates are 0 (existing
  behaviour) and the Stan model.

## Results

Bit-identical. Checked with `dev/compare-commits.R` (installed devel
against the branch), which covers the `decontPro` and `decontProAmbient`
Stan paths as well as the decontX paths.

## Tests

A pure refactor, so no new test can fail first. The hand-built-fit oracle
test in `test-decontPro.R` pins all three matrices and the parameter list,
and the real-fit test checks that the three parts add up to the counts;
both must pass unchanged.

## Benchmark

`decontPro()` is dominated by Stan, so the whole run won't show the
change. `.process_stan_vb_out()` alone is timed with `bench::mark` (time
and R-allocated memory) on a synthetic hand-built fit at 40 ADTs x 50,000
droplets and at 215 x 52,222 (the size in issue #46).

## Other

Version 1.11.22. One NEWS line only if the memory saving is worth
mentioning. ROADMAP I-e notes it as done. No doc or vignette changes;
decontPro isn't reached through celda or singleCellTK.

## Measurements

`bench::mark`, 10 iterations, devel's and the branch's
`.process_stan_vb_out()` sourced into one R 4.6.1 session (macOS arm64),
on a synthetic hand-built fit (8 cell types, Poisson counts). Outputs
`identical()` at both sizes.

| ADTs x droplets | Median, devel | Median, branch | R alloc, devel | R alloc, branch |
|---|---|---|---|---|
| 40 x 50,000 | 49 ms | 40 ms | 321 MB | 305 MB |
| 215 x 52,222 | 267 ms | 210 ms | 1.76 GB | 1.67 GB |

The allocation saving is one N x M double matrix at both sizes (16 MB and
90 MB). The time saving is small next to the Stan fit, so there is no
NEWS line.
