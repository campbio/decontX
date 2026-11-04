# ADR 0004: Separate a container-agnostic core from the decontX and decontPro methods

## Status

Accepted (2026-11-04). Proposed as a stub on 2026-10-27. The decisions
below were made during the Phase 0 review (2026-10-18 to 2026-10-27);
the points that stub left open were decided on 2026-11-03 (see "Decided
2026-11-03").

## Context

The deep review (`dev/audits/2026-10-deep-review.md`, area 4) found that
the core of decontX is tangled with the SingleCellExperiment method:

- The decontX argument surface is repeated at about 11 sites, with real
  default drift (`maxIter` 500 vs 200, `varGenes` 5000 vs NULL vs 2000).
  Internal calls pass every argument, so no user sees the drift today
  (A4-01, D1-20).
- Validation is scattered and repeated per batch, which is the root of
  the `z` length bug (A4-02, A1-01).
- `.decontX` is 307 lines and `.decontXoneBatch` 190 (A4-03, D1-21).
- `runParams` is used as data, stores the processed batch, and lacks
  `seed` and `iterLogLik` (A4-04, A1-26).
- Core logic lives in the SCE method: UMAP assembly and background
  deduplication (A4-05). There are three copies of the seed block
  (A4-06), and `z` is duplicated inside `estimates` (A7-22).
- Batched runs assemble the result with about 7.5 full copies of the
  counts (A5-05).
- HDF5 input spends 63% of its runtime writing the result back to a
  temporary HDF5 file, and the saved object then points at that file
  (A5-06).

Phase II adds Seurat, AnnData and ExperimentSubset methods
(ADR-0005). Each would otherwise repeat the logic now in the SCE
method. A core that takes matrices and returns plain lists lets every
method stay a thin read-call-write wrapper.

## Decision

We will split decontX and decontPro into a container-agnostic core and
thin per-class methods, keeping results bit-identical.

**Core contract.** `.decontX(counts, z, batch, background, bgBatch,
params)` takes counts (genes x cells, dense or sparse), `z`, `batch`,
background counts, background batch, and the params list, with no
defaults. It checks the data and handles the background (shape checks,
removing barcodes that are also in the counts, matching rows to genes).
It returns a plain list: `runParams`, `estimates`, `decontXcounts`,
`contamination` and `z`, as today, then `batch`, the per-cell batch
labels used for the run. The core never reads or writes a container.
Decided 2026-10-18 (area 4, decision 4): the `ANY` method returns the
extra fields, with a NEWS entry.

**`runParams` is a record, not data.** It keeps today's elements and
values: the raw `z`, the processed per-cell `batch` and
`batchBackground`, and the parameters as given, including `seed` and
`iterLogLik`. `batchBackground` is always present, `NULL` without a
background (#96). Methods read `result$batch`, never `runParams`.

**UMAPs.** Each batch's UMAP stays in `estimates[[b]]$UMAP`, with that
batch's cells only; the core returns no separate UMAP element. A shared
helper, `.decontXUMAPs()`, returns the matrices for containers that need
one row per cell (SCE `reducedDims`, AnnData `obsm`): the single batch's
UMAP as it is, or each batch's padded with `NA` for the other batches'
cells. Seurat's per-batch reductions take `estimates[[b]]$UMAP`
directly.

**decontPro.** `.decontPro()` is already decontPro's core (matrices in,
a list out) and keeps returning today's list (`decontaminated_counts`,
`background_counts`, `ambient_counts`, `parameters`). ADR-0005 decides
what the container methods write and return, and ADR-0006 the argument
names. A decontPro params constructor waits for ADR-0006.

**One params constructor.** `.decontXParams()` validates and builds the
parameters, with no defaults of its own, and returns them as a plain
list. The user-facing methods keep full signatures with their defaults,
and a test checks that the methods' formals agree with each other and
with `.decontXParams()` (pattern A, decided 2026-10-18, A4-01). New
decontX arguments are added to `.decontXParams()` and to every method
signature, and the formals test covers them; they don't need an ADR of
their own.

**Decomposition.** `.decontX()` and `.decontXoneBatch()` are split into
named helpers in the same file (A4-03). The core does the background
deduplication and the seeding once, and the UMAP padding leaves the SCE
method for `.decontXUMAPs()` (A4-05, A4-06). `z` is no longer duplicated
inside `estimates` (A7-22).

**Batched assembly.** One `cbind` and one reorder replace the
per-batch copies, which cuts peak copies from about 7.5 to about 2 with
identical results (A5-05).

**HDF5 input.** HDF5 input returns an in-memory sparse result wrapped as
a DelayedMatrix, instead of writing back to a temporary HDF5 file. The
class stays DelayedMatrix and the values are identical; users save to
disk explicitly (for example `saveHDF5SummarizedExperiment()`).
HDF5Array is not added to Suggests. Decided 2026-10-19 (A5-06), which
reverses the earlier "HDF5 in, HDF5 out" decision. NEWS.

**File moves last.** After the refactor, a pure file-move PR moves
`decontXcounts` and `simulateContamination` into their own files; the
core stays in `R/decon.R` (A4-17, decided 2026-10-18).

**Bit-identical.** Committed reference outputs (tolerance 1e-10 on the
user-`z` paths) and the `dev/` cross-commit `identical()` comparison
(no-`z`, `legacyInit`, decontPro) from Phase I-a must pass unchanged on
the single-batch, batched, background and user-`z` paths.

Alternative considered: pattern B, passing `...` through the internal
calls instead of a constructor. It lost because it hides the argument
surface and keeps the default drift possible.

## Consequences

- Phase II methods (ADR-0005) become thin: read the inputs, call the
  core, write the results back. Logic in a method signals a gap in the
  core.
- One place validates parameters, which closes the class of bugs behind
  A1-01 and A1-10.
- The `ANY` method's returned list gains `batch`, and `runParams`
  always has `batchBackground`. Code that checks the list's length or
  names may notice. celda forwards to `decontX::decontX()` and reads
  nothing back; singleCellTK reads `runParams$batch` and
  `$batchBackground`, which don't change. The SCE results are identical
  apart from `runParams$batchBackground`.
- HDF5 users get an in-memory result, so a very large HDF5 input now
  needs memory for the sparse result. NEWS says so.
- Batched runs use less memory.
- The refactor lands as several PRs (core and params, decomposition,
  batched assembly, HDF5, file move), each held to the reference tests.

## Decided 2026-11-03

The stub's open points, and the order of the PRs:

- **New arguments from Phase I-e** (the initialization opt-out when `z`
  is given, A5-01, and `numThreads`, A5-02): recorded in #79's spec, not
  in this ADR or one of their own. This ADR only states the rule for new
  arguments (see "One params constructor").
- **List element names:** as in "Core contract". The UMAPs stay in
  `estimates` (see "UMAPs"). A padded `UMAP` element in the result was
  rejected because it duplicates every UMAP in the matrix result.
- **decontPro:** no change in #76 (see "decontPro").
- **`runParams`:** keeps its contents (see "`runParams` is a record").
  Storing the raw `batch` was rejected: it changes saved metadata and
  the per-cell vectors singleCellTK splits by sample
  (`R/celda_decontX.R` there).
- **PRs, each from and into devel, "Part of #76":** (1) the core
  contract, `.decontXParams()`, the formals test, the logic moved out of
  the SCE method, and #96; (2) decomposition, one seed block, dead code,
  `.logMessages`, `z` out of `estimates`; (3) batched assembly; (4) HDF5;
  (5) the file move, which closes #76. No `refactor/` base branch: each
  piece gives the same results and leaves devel releasable.
