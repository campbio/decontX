# Less memory for decontPro's split of the counts (#111)

Approved in chat 2026-11-10 (variant E of the evaluation). Follows A4-20
(#80 part 3, `dev/plans/2026-11-08-stan-process-tidy.md`).

## Problem

`.processStanVbOut()` in `R/stan_helpers.R` splits decontPro's counts
into native, ambient and background parts. It builds two N x M double
matrices (ADTs x droplets) only to scale by a per-droplet and a per-ADT
vector (`scaling_factor` from `OC`, `ambient_profile` from `p`), and
keeps every intermediate alive until the three splits are done. At the
size in #46 (215 x 52,222) each such matrix is about 86 MB, and the
function's peak sets decontPro's memory after the Stan fit.

## Change

- Leave `OC` out of the rates. Each of a droplet's three rates is
  `OC[m]` times a term, so `OC[m]` cancels in every rate / total share.
  The scaling matrix goes.
- `dat$p` recycles down each droplet's column instead of an N x M
  profile matrix.
- `1 - background_est` is computed once; the background rate is
  `background_est` itself.
- One `counts / total` matrix is shared by the three splits
  (`counts_per_rate * rate`, so the counts' dimnames lead), and the
  profile and rate matrices are removed with `rm()` once used.
- `drop = FALSE` on `r_est[dat$cell_type, ]` keeps the native profile a
  matrix for any shape (the second bug noted in #110; unreachable once
  #110 requires 2 droplets).
- Outputs, their names and the parameter list are unchanged. Out of
  scope: the Stan model and the 0/0 case when all three rates are 0.

## Results

Not bit-identical, as agreed: removing `OC` and sharing `counts / total`
changes the order of the floating-point operations. The parameter list is
identical. On the synthetic fits below the three splits move by at most
7.8e-16 relative, and the three parts add up to the counts as closely as
before (max |sum - counts| 3.6e-15 on both). `dev/compare-commits.R`
(devel against the branch) reports the decontX paths identical and the
decontPro paths within the same tolerance. The committed references
(`decontx-refs.rds`) cover decontX only, so none change.

## Tests

The hand-built-fit oracle test in `test-decontPro.R` (exact splits for 2
ADTs x 3 droplets) and the real-fit test that the parts add up to the
counts must pass unchanged. No new test: this is a refactor with no new
behaviour, and the oracle test fails on any wrong formula.

## Benchmark

Devel's and the branch's `.processStanVbOut()` sourced into one R 4.6.1
session (macOS arm64), on a synthetic hand-built fit (8 cell types,
Poisson counts). `bench::mark`, 10 iterations, for time and R-allocated
memory; peak RSS from `/usr/bin/time -l` in a fresh process per version,
minus a run that builds the fit and calls nothing (3 runs each).

| ADTs x droplets | Median, devel | Median, branch | R alloc, devel | R alloc, branch |
|---|---|---|---|---|
| 40 x 50,000 | 69 ms | 24 ms | 306 MB | 199 MB |
| 215 x 52,222 | 299 ms | 162 ms | 1.67 GB | 1.09 GB |

Devel's median moves between runs (35 to 78 ms at 40 x 50,000 across
this session, 228 to 299 ms at 215 x 52,222); the branch was faster in
every run, by about 30% to 65%.

Peak RSS at 215 x 52,222: baseline 623 to 637 MB; devel 1,569 to 1,576
MB (about 945 MB above baseline); branch 1,133 to 1,146 MB (about 510 MB
above baseline). The split's own peak drops by about 46%, from about 11
to about 6 of the 86 MB matrices.

`dev/compare-commits.R campbio/devel HEAD`, run with the split commit
at HEAD: all 13 decontX paths identical; `decontPro` and
`decontProAmbient` differ only in the three splits, with mean relative
differences of 1.6e-16 to 2.1e-16.

## Other

Version 1.11.29. A NEWS line, since decontPro's memory after the fit
drops and results can differ in the last digits. ROADMAP I-e marks #111
done. No doc or vignette changes; decontPro isn't reached through celda
or singleCellTK.
