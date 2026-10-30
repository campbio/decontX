# ADR 0004: Separate a container-agnostic core from the decontX and decontPro methods

## Status

Proposed (2026-10-27). A stub for maintainer review: the decisions below
were made during the Phase 0 review (2026-10-18 to 2026-10-27); the open
points at the end are not yet decided.

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

**Core contract.** `.decontX()` takes counts (genes x cells, dense or
sparse), `z`, `batch`, background counts, background batch, and a params
object. It returns a plain list: decontaminated counts, contamination,
cluster labels, batch, estimates, the UMAP, and `runParams`. `runParams`
holds the validated parameters plus `seed` and `iterLogLik`. The raw `z`
is kept. `.decontPro()` follows the same pattern. The core never reads
or writes a container. Decided 2026-10-18 (area 4, decision 4): the
`ANY` method returns the extra fields, with a NEWS entry.

**One params constructor.** `.decontXParams()` validates and builds the
parameters, with no defaults of its own. The user-facing methods keep
full signatures with their defaults, and a test checks that the
methods' formals agree (pattern A, decided 2026-10-18, A4-01).

**Decomposition.** `.decontX()` and `.decontXoneBatch()` are split into
named helpers in the same file (A4-03). The core does the UMAP assembly,
background deduplication and seeding once (A4-05, A4-06). `z` is no
longer duplicated inside `estimates` (A7-22).

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
- The `ANY` method's returned list gains fields (batch, UMAP, a fuller
  `runParams`). Code that checks the list's length or names may notice.
  celda forwards to `decontX::decontX()` on SCE and is unaffected by the
  SCE path, which stays bit-identical.
- HDF5 users get an in-memory result, so a very large HDF5 input now
  needs memory for the sparse result. NEWS says so.
- Batched runs use less memory.
- The refactor lands as several PRs (core and params, decomposition,
  batched assembly, HDF5, file move), each held to the reference tests.

## Open points

- Whether the two new arguments from Phase I-e, the opt-out of
  initialization when `z` is given (A5-01) and the opt-in `numThreads`
  (A5-02, default 1), are recorded here or in a short ADR of their own.
  Their behaviour is decided (2026-10-19): initialization still runs by
  default, and `numThreads > 1` changes PCA by about 1e-11. Their names
  are not.
- The exact names of the list elements returned by `.decontX()` and
  `.decontPro()`.
- How `.decontPro()`'s list relates to the current `decontPro()` matrix
  return, which ADR-0006 also changes.
