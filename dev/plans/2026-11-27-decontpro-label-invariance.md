# decontPro results independent of the cell-type label names (#148)

Design approved in chat on 2026-11-27, with both parts. The `r` start
was revised to uniform during implementation (approved 2026-11-27).
Version 1.11.41, branch `fix/decontpro-label-invariance` from `campbio/devel` 6abc6b9.

## Problem

`.recodeCellType()` codes the labels 1..K in sorted order, or in factor
level order. `decontPro()` gives starting values only for `delta` and
`background`, so Stan draws random starts for `r` (one native profile per
cell type, in code order), `delta_mean`, `background_mean`, `tau_a` and
`tau_b`. Renaming a label can reorder the codes, so each cell type gets
a different random start. With coarse labels, ADVI then settles in a
different optimum: in #148, renaming "mixed" to "unclear" switched the
fit from removing CD8a to removing CD4.

Fixed starting values alone don't fix this. At every gradient step, ADVI
draws a standard-normal vector over all unconstrained parameters, and
`r`'s K blocks sit in that vector in code order. A reorder gives each
cell type different noise, even from the same start.

## Decisions

1. **Stan codes by first appearance.** `.decontPro()` codes the cell
   types in the order they first appear in the droplets (the columns of
   `filtered_counts`), which doesn't depend on the label names or on a
   factor's level order. Renaming the labels or reordering the levels
   gives Stan identical input, so the results are identical. The rows of
   `parameters$r_est` are put back in the order they have now (sorted
   labels, or a factor's levels with unused levels dropped) and named as
   now. The `cellType` documentation stays as it is.
2. **Fixed starting values for every parameter**, so Stan draws no
   random start:
   - `r[k]`: the uniform profile, 1/N for each ADT (the centre of the
     simplex; revised during implementation, see below);
   - `delta` = 1e-4 and `background` = 1e-2, as now; `delta_mean` = 1e-4
     and `background_mean` = 1e-2, matching them;
   - `tau_a` = `tau_b` = 1.

   They come from one helper, `.decontProInits(N, M, K)`.
3. **No option to keep the old behaviour** (maintainer, 2026-11-27).
   Results change for every user; NEWS says so, and that decontX 1.11.40
   or earlier gives the old numbers. The old results depended on an
   arbitrary label order, and decontPro already broke compatibility in
   this cycle (ADR-0006).
4. **No ADR.** The Stan model (`inst/stan/shrinkage.stan`) is unchanged;
   only the data order and the starting values change.

## Probe behind the decisions

The PBMC article's data (8,000+ cells, 17 ADTs). The fine labels are the
article's 8 types; the coarse labels have one "T" label. "current" is
decontX 1.11.40; "chosen" is decisions 1 and 2. Each value is the share
of an ADT's counts removed.

| Fit | Labels | Seed | CD4 | CD8a | CD45RA | Ambient share | ELBO |
|---|---|---|---|---|---|---|---|
| current | fine | 12345 / 1 | 0.20 | 0.19 | 0.32 / 0.31 | 0.189 | -5.154M / -5.155M |
| chosen | fine | 12345 / 2 | 0.20 | 0.19 | 0.32 | 0.189 | -5.157M / -5.155M |
| current | coarse "mixed" | 12345 | 0.28 | 1.00 | 0.55 | 0.337 | -10.68M |
| current | coarse "unclear" | 12345 | 0.96 | 0.18 | 0.46 | 0.372 | -10.60M |
| current | coarse | 1, 2 | 0.96-0.97 | 1.00 | 0.47-0.56 | 0.42 | |
| chosen | coarse, both names | 12345 | 0.28 | 1.00 | 0.55 | 0.337 | -10.68M |
| chosen | coarse | 2 | 0.28 | 1.00 | 0.55 | 0.336 | -10.69M |
| chosen | coarse | 1 | 0.97 | 1.00 | 0.56 | 0.424 | -10.66M |

- With the chosen start, "mixed" and "unclear" give `identical()`
  results.
- With the fine labels, the chosen start fits as well as the current
  code and the results barely move.
- With the coarse labels, the result still depends on the seed: ADVI's
  noisy gradient steps on a misspecified model reach optima with similar
  ELBOs. No start fixes that. The fine-label fits are stable across
  seeds.

The decontPro vignette's 500 cord blood droplets (`set.seed(42)`), share
of each ADT removed:

| Fit | Ambient | Background | CD8 | CD16 | CD56 | Correlation with the all-droplet fit |
|---|---|---|---|---|---|---|
| current, default priors | 0.084 | 0.012 | 0.08 | 0.18 | 0.19 | |
| chosen, default priors | 0.102 | 0.012 | 0.09 | 0.20 | 0.20 | |
| current, 5x priors | 0.196 | 0.042 | 0.26 | 0.48 | 0.49 | 0.989 |
| chosen, 5x priors | 0.215 | 0.042 | 0.27 | 0.49 | 0.52 | 0.983 |
| current / chosen, 10x priors | 0.390 / 0.390 | 0.057 / 0.059 | 0.37 / 0.37 | 0.62 / 0.62 | 0.59 / 0.60 | |

The correlation is between the log counts decontaminated with 5x priors
and those of a default fit on all 7,700 labelled droplets; the vignette
now says 0.98.

Rejected:

- **`r` at each cell type's observed ADT proportions** (+1 per ADT), the
  start first approved. Those proportions include the contamination, so
  the native profile starts out explaining nearly every count. The fits
  settled in a worse optimum: fine labels ELBO -5.35M against -5.15M,
  CD45RA removed 0.22 against 0.32; on cord blood, 54% of CD16 removed
  with the default priors against 18%, and a correlation of 0.968. Its
  coarse fits removed CD8a with all 3 seeds.
- **`tau_a` and `tau_b` at the prior scale** (`deltaSd`, `backgroundSd`).
  It pins every droplet's `delta` to its mean; ADVI stopped early with
  much worse ELBOs (fine labels -5.34M; coarse labels down to -16.9M,
  one fit with no ambient contamination). Its coarse fits kept CD4 and
  CD8a, which looked right, but the model itself rates them worse.

Pitfall found on the way: rstan silently ignores an `r` start given as a
list of vectors, and starts `r` at random. The first probe passed it that
way, so its "fixed start" fits had a random `r`. `r` must be a K x N
matrix; a unit test pins the form and a real fit with one cell type
checks that rstan accepts it.

## Changes

- `R/decontPro.R`: `.recodeCellType()` returns the first-appearance
  codes, the label of each code, and the output order; `.decontPro()`
  passes all starting values from `.decontProInits()`.
- `R/stan_helpers.R`: `.processStanVbOut()` computes the splits with the
  Stan-order `r` and stores `r_est` in the output order.
- `?decontPro`, `seed`: when the labels don't fit the data well (for
  example one label for two cell types whose markers differ), the
  result can depend on the seed, so running a second seed is a way to
  check.
- `vignettes/articles/decontPro_pbmc10k.Rmd.orig`: one sentence in the
  coarse-label section saying which marker is lost depends on the seed.
  The `lost_marker` check also requires the other marker to stay below
  0.8, so the text's "one kind of T cell" stays true. The maintainer
  re-knits it (`make knit-article FILTER=decontPro_pbmc10k`); with seed
  12345 the coarse fit still loses CD8a.
- NEWS (1.11.41), DESCRIPTION version.
- `vignettes/decontPro.Rmd`: the correlation with the all-droplet fit
  becomes 0.98. Its statements about the density plots rest on shares
  that move by at most 0.03; checked in the rendered vignette.

## Tests (`tests/testthat/test-decontPro.R`, 10 x 8 fixture)

Written first; each fails on its assertion with the current code.

- Renaming the labels (with a rename that changes their sorted order)
  gives `identical()` decontaminated, ambient and background counts,
  and `r_est` equal up to the row names.
- A factor with reversed levels gives identical counts, with the `r_est`
  rows in the reversed level order.
- `.decontProInits()`: `r` is a K x N matrix of 1/N; the other
  parameters have the right shapes and values.
- `decontPro()` passes `.decontProInits()` to Stan and codes the cell
  types by first appearance (mocked Stan call).
- A real fit with one cell type (one-row `r`).
- The existing `r_est` row-name tests keep passing unchanged.
- No decontPro references are in `decontx-refs.rds`, so none need
  regenerating; `dev/compare-commits.R` shows which decontPro paths move.

## Follow-up (not in this PR)

An issue proposing several seeds per fit and keeping the fit with the
best ELBO, to address the coarse-label instability itself.
