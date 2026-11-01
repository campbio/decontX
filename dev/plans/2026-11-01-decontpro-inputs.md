# decontPro: split Seurat layers, `cell_type` labels, unknown arguments (#72)

Status: approved by the maintainer, 2026-11-01 (recode `cell_type` and
name the `r_est` rows). Roadmap Phase I-b; audit findings A1-04, A1-11
and A1-13. A1-05 (reading ADT rather than RNA) is left for the Phase II
container methods.

## Why

Reproduced on devel at f380a36 (1.11.8) with the 10 x 8 fixture.

- **Split layers (A1-04).** `LayerData(x, layer = "counts")`
  (`decontPro.R:97, 99`) returns the first of several counts layers with
  a warning. A Seurat v5 object split in two passes 10 x 4 to
  `.decontPro`: a full-length `cell_type` then fails in Stan, and a
  half-length one runs on half the cells.
- **`cell_type` and counts (A1-11).** Nothing is checked. Gaps (1, 3, 5,
  7) give "no more scalars to read"; NA gives "argument is not a matrix";
  a wrong length, character, factor, 0-based or non-integer labels, and
  non-integer counts make `rstan::vb()` return a list, so
  `.process_stan_vb_out()` stops with "no applicable method for `@`
  applied to an object of class list". Gaps come up in practice:
  `as.integer(Idents(seu))` after subsetting.
- **`...` (A1-13).** Unknown arguments are dropped: `seed = 99` gives the
  default result, a `delta_SD = 1e-2` typo runs with the default
  `delta_sd`, and `batch` is ignored.

## What changes

1. **Seurat counts.** A helper reads the counts of a Seurat object's
   default assay. A v5 assay with several counts layers is joined with
   `SeuratObject::JoinLayers()` on the extracted copy (the user's object
   is untouched), so columns follow `Cells()`. Used for `filtered_counts`
   and `ambient_counts`.
2. **`cell_type`.** One entry per droplet and no NA, else a stop naming
   the argument. Any labels are accepted and recoded to 1..K for Stan:
   sorted unique values, or a factor's levels with unused levels
   dropped. Valid 1..K input gets the same codes, so its numbers don't
   change. `parameters$r_est` gets the labels as row names.
3. **`filtered_counts`.** Must be numeric whole numbers, with no NA and
   no negatives (the model declares `int<lower = 0>`).
4. **Stan failure.** If `.call_stan_vb()` doesn't return a `stanfit`,
   decontPro stops with "Stan's variational inference failed" instead of
   the `@` error.
5. **`...`.** Any extra argument, named or not, stops in all three
   methods, naming the arguments.
6. **Docs.** `@param cell_type` (labels and the `r_est` row order),
   `@param filtered_counts` (whole numbers; split layers are joined),
   `@param ...` (unused; extra arguments are an error), `@return`
   (`r_est` row names).

## Tests

Written first; each fails on its assertion on devel.

- A Seurat object split into two layers with interleaved cells (and a
  split `ambient_counts` object) passes all 8 columns, in `Cells()`
  order, to a mocked `.decontPro`, with no warning.
- NA, short and long `cell_type` stop with the argument named.
- Gapped, character, factor (with an unused level) and 0-based labels
  reach a mocked `.call_stan_vb` as codes 1..K with the right K.
- A real 10 x 8 fit with gapped labels equals the 1..K fit, apart from
  the `r_est` row names, which are the labels.
- Non-integer, negative and NA counts stop with clear messages.
- A mocked `.call_stan_vb` returning a list stops with the new message.
- `seed`, `delta_SD`, `batch` and an unnamed extra stop in the matrix,
  SingleCellExperiment and Seurat methods.

## NEWS, version, release

- NEWS entries for the three fixes and the `r_est` row names.
- Version 1.11.10 (PR #92 takes 1.11.9); rebased locally on devel once
  #92 merges.
- Devel only: crashes with unclear messages; the split-layer half-run
  warned.

## Effect on results and related packages

decontPro numbers don't change for input that ran before; `r_est`
gains row names, so `dev/compare-commits.R` reports the two decontPro
paths as differing in that attribute only. celda and singleCellTK
don't call decontPro.
