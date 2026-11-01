# Fix crashes on small batches and some input classes (#71)

Status: approved by the maintainer, 2026-11-01. Roadmap Phase I-b; audit
findings A1-06, A1-08, A1-14, A1-15 and D2-02.

## Why

Reproduced on devel at f380a36 (1.11.8).

- **Small batches crash (A1-08).** `counts[, batch == bat]`
  (`decon.R:592`) drops a one-cell batch to a vector, and
  `estRmat[, !(allCol %in% subCol)]` (`:682`) does the same when the
  other batches hold one cell in total: "subscript out of bounds", or
  "no method or default for coercing 'numeric' to 'dgCMatrix'" for a
  one-cell `x`. `drop = FALSE` alone isn't enough, because the
  initialization always runs PCA and UMAP:
  - default initialization: 3 cells run, 2 cells stop in the PCA
    ("requested number of singular values cannot be greater than the
    smaller matrix dimension");
  - `legacyInit = TRUE`: 15 cells run, 14 stop in uwot ("n_neighbors
    must be smaller than the dataset size").

  A one-cell batch with `z` and no `background` already stops cleanly
  with the one-cluster error from #69.
- **DelayedMatrix with a dense seed (A1-06).** The seed's class is
  `c("matrix", "array")`, so `seed.class == "HDF5ArraySeed"` (`:737`)
  is a length-2 condition: "the condition has length > 1", after the
  whole EM has run.
- **Input classes (A1-15).** `as(x, "dgCMatrix")` (`:604`, `:608`) is
  deprecated in Matrix: dgRMatrix stops ("no method or default for
  coercing"), dgTMatrix warns. A data.frame stops with the same coercion
  error; a SummarizedExperiment `x`, or a matrix `x` with a
  SingleCellExperiment `background`, stops with "invalid 'type' (S4) of
  argument"; a SingleCellExperiment `x` with a matrix `background` stops
  with "unable to find an inherited method for function 'assay'". The
  only class check is `is.null(dim(x))`, and its message ("At least 2
  genes need to have non-zero expressions") doesn't describe the problem.
- **The output class (A1-14).** The `try()` block at `:749-757` meant to
  give back the input's class. For dense input `canCoerce()` gets a
  length-2 class and errors inside the silent `try()`, so dense input
  returns a dgCMatrix (kept, per the maintainer's 2026-10-18 decision).
  It does still convert other classes back: dgTMatrix input returns a
  dgTMatrix and SVT_SparseMatrix input returns an SVT_SparseMatrix.
- **C++ guard (D2-02).** `decontXInitialize()` (`src/DecontX.cpp:227`)
  is the one EM entry point without a `min(z) < 1` check: 0 corrupts
  memory, negatives segfault. R code never passes such a `z` today.

## What changes

1. **Minimum batch size.** Before any batch runs, `.decontX()` stops if
   a batch has fewer than 3 cells (15 with `legacyInit = TRUE`). The
   message names the batch and its size, or says `x` when there's no
   `batch`. Both subsets get `drop = FALSE` as well.
2. **Input classes.** `x` (the ANY method, and the SCE method's assay)
   and `background` must be a numeric base matrix, a Matrix-package
   matrix, a DelayedMatrix, or a SparseArray matrix (such as
   SVT_SparseMatrix). Anything else stops before any work, naming the
   argument and its class; a SummarizedExperiment `x` is told to convert
   with `as(x, "SingleCellExperiment")`. With a SingleCellExperiment `x`,
   `background` must be a SingleCellExperiment too, as the docs already
   say ("the same data format as x").
3. **One conversion.** A helper converts each batch to a dgCMatrix with
   `as(x, "CsparseMatrix")`, then `generalMatrix` and `dMatrix`, so
   dgR/dgT/SVT/DelayedMatrix/symmetric/logical input converts without
   errors or deprecation warnings. Already-dgCMatrix input is untouched.
4. **Output class.** The `try()` coercion is deleted: `decontXcounts` is
   a dgCMatrix, or a DelayedMatrix for DelayedArray input. An
   HDF5ArraySeed still gives an HDF5Matrix; every other seed gives
   `DelayedArray(<dgCMatrix>)`. The seed check uses
   `methods::is(seed, "HDF5ArraySeed")` over `DelayedArray::seedApply()`.
5. **C++.** `decontXInitialize()` gets `min(z) < 1` with the message of
   the other entry points ("The entries in 'z' need to be ... "). NA
   (INT_MIN) is caught by the same check.
6. **Docs.** `@param x`/`background` list the accepted classes and the
   minimum batch size; `@return` says `decontXcounts` is a dgCMatrix, or
   a DelayedMatrix for DelayedArray input.

## Tests

Written first; each fails on its assertion on devel.

- A one- and a two-cell batch, and a two-cell `x`, stop with the
  minimum-size message naming the batch; three cells run. With
  `legacyInit = TRUE` (skipped without scater), 14 cells stop and 15 run.
  A run whose other batches total one cell no longer hits "subscript out
  of bounds" (the size check now catches it; `drop = FALSE` is covered
  by a three-cell batch next to a large one).
- dgRMatrix, dgTMatrix and SVT_SparseMatrix input run without warnings
  and give the dense-input contamination and a dgCMatrix
  `decontXcounts` equal to the dense-input one.
- A DelayedMatrix over a dense seed returns a DelayedMatrix whose values
  equal the dense-input result.
- data.frame and SummarizedExperiment `x`, a SingleCellExperiment
  `background` with a matrix `x`, and a matrix `background` with a
  SingleCellExperiment `x` each stop with a message naming the argument
  and class.
- `decontXInitialize()` called directly with a 0, a negative and an NA
  in `z` stops with the range message.
- The comment in the existing "sparse and DelayedArray input" test that
  says "#71 changes this" is updated: #71 keeps dense in, dgCMatrix out.
- No change elsewhere: the stored references still match, and
  `dev/compare-commits.R devel HEAD` reports every path identical.

## NEWS, version, release

- NEWS: small batches stop with a clear message; DelayedMatrix with a
  dense seed, dgRMatrix and dgTMatrix input work; unsupported classes
  stop with a clear message; `decontXcounts` is always a dgCMatrix (or
  DelayedMatrix), so dgTMatrix and SVT_SparseMatrix input no longer come
  back in their own class.
- Version 1.11.9.
- Devel only: crashes on uncommon inputs, no wrong results.

## Effect on related packages

celda forwards to `decontX()`. Inputs that crashed now run or stop
clearly; dgTMatrix and SVT_SparseMatrix input return a dgCMatrix. No
numbers change for any input that ran before.
