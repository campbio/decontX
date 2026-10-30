# ADR 0005: Add Seurat, AnnData and ExperimentSubset methods that return the same class

## Status

Proposed (2026-10-27). A stub for maintainer review: the decisions below
were made during the Phase 0 review (2026-10-18 to 2026-10-27) and the
area 6 spike; the open points at the end are not yet decided. Depends
on ADR-0004 (the container-agnostic core).

## Context

decontX and decontPro dispatch on `SingleCellExperiment`, `Seurat`
(decontPro only) and `ANY`. Users of other containers convert by hand,
and the vignette's Seurat round trip is broken (A7-02). The review found
more:

- decontPro on a Seurat object with split layers reads only the first
  layer (A1-04), and it reads the default assay or the SCE `counts`
  assay, so it ran on RNA instead of ADT; `altExpName` is swallowed by
  `...` (A1-05, A6-09, A4-12).
- The plots read containers directly, so they can't accept other
  classes, and the marker plots subset the whole SCE (A4-18, A4-19). The
  contamination plot depends on UMAP column names, which an h5ad round
  trip drops (A6-04).
- `setOldClass("AbstractAnnData")` alone doesn't dispatch S4 methods for
  anndataR objects (A6-03).
- The three ExperimentSubset classes share no class or predicate, so
  `is(x, "SubsetSingleCellExperiment")` misses the spatial and tree
  classes (A6-02). ExperimentSubset also had assay-name scoping and
  colData bugs.

The area 6 spike showed that every target class can be supported
without a temporary SCE object, using the core from ADR-0004.

## Decision

We will add one method per class on `decontX()`, `decontPro()` and
`decontXcounts()`/`decontXcounts<-`. Each reads its inputs from the
object, calls the core, and writes the results back into an object of
the same class.

**Conventions for every class**
- `z` and `batch` accept a vector or a single string naming a column in
  the object's cell metadata (`colData`, `meta.data`, `obs`, or the
  subset's colData).
- Accessors return genes x cells for every class.
- decontPro container methods return the updated object: a clean break,
  marked Breaking in NEWS. The `ANY` method keeps returning the list.
  The decontPro break and the new methods land in the same devel cycle,
  so one Bioc release carries the whole API change.
- The plots take any supported object through an internal results
  getter that extracts matrices first, so they no longer copy every
  assay (A4-18, A4-19). They index UMAP columns by position (A6-04).

**Read and write locations** (verified in the spike; audit area 6):

| Class | Input | Decontaminated counts | Per cell | UMAP | Run info |
|---|---|---|---|---|---|
| SCE, SpatialExperiment | `assay(x, assayName)`, validated against `assayNames(x)`; decontPro `assay(altExp(x, altExpName), assayName)` | assay `decontXcounts`; decontPro assay `decontPro` in the altExp or main | colData | `decontX_UMAP`, or `decontX_<b>_UMAP` per batch, padded with NA | `metadata(x)$decontX`; decontPro in the altExp or main metadata |
| ExperimentSubset | `assay(x, subsetName)` | `setSubsetAssay(x, s, m, "decontXcounts")` | `colData(x, subsetName = s) <-` | `reducedDim(x, , subsetName = s) <-` | `metadata(x, subsetName = s)$decontX` |
| Seurat | Assay5: `JoinLayers` on a copy, then `LayerData` counts; v3: `LayerData` | new assay `decontXcounts` | `meta.data` by `Cells()` | reduction `decontX_UMAP`, key `decontXUMAP_`; per-batch reductions hold only that batch's cells | `Misc(seu, "decontX")` |
| AnnData (concrete, not views) | `t(X)` or `t(layers[[l]])` | `layers[["decontXcounts"]]` | `obs` | `obsm[["X_decontX_umap"]]`, per batch `X_decontX_<b>_umap` | `uns$decontX` (run parameters and per-batch estimates only) |

**Seurat.** Split count layers are joined on a copy, in `Cells()` order.
A supplied `batch` is used as-is; if `batch` is NULL and the layers are
split, the layer split becomes the batch, with a message. Results go to
a separate `decontXcounts` assay (decided 2026-10-19/20). v3 assays have
no `JoinLayers`; reduction keys are sanitized; `Misc` is cleared before
a re-run (A6-08). decontPro reads the ADT assay.

**AnnData.** A thin wrapper over anndataR (Suggests). It extracts the
matrix, resolves `z`/`batch`, calls the core (which clusters when `z` is
NULL), and writes back to `layers`/`obs`/`obsm`/`uns`, with no temporary
SCE. AnnData is R6, so the object is modified in place; the docs show
`adata <- decontX(adata)` and say the original changes too. Dispatch is
an `inherits(x, "AbstractAnnData")` redirect in the `ANY` methods, not
`setOldClass` (amendment, 2026-10-27). anndataR's own error for views
is let through; a read-only backed object gets a clear error with a
`mode = "r+"` hint where the mode can be queried (A6-10). The
h5ad -> SCE -> h5ad route through anndataR is documented as an
alternative. zellkonverter was rejected: it needs a conda Python
environment at runtime and returns an SCE, not an AnnData object.

**ExperimentSubset: the direct approach** (decided 2026-10-27, A6-12 to
A6-19). Read `assay(x, subsetName)` and the subset's colData (with
`parentColData = TRUE` for `z`/`batch` lookup), call the core, and write
back with `setSubsetAssay(..., unlink =)`, `colData(x, subsetName =) <-`
(the subset's own columns only), `reducedDim(x, name, subsetName =) <-`
and `metadata(x, subsetName =) <-`. Not `getSubsetExperiment()` /
`setSubsetExperiment()`.
- `unlink = TRUE` is needed only when the main object already has
  `decontXcounts` (A6-14).
- Writing back only the subset's own columns avoids the silent drop of
  overwritten inherited columns (A6-13).
- The subset's primary assay is the default input; `assayName =
  "decontXcounts"` is rejected (A6-15).
- `colData(..., parentColData = TRUE)` can return duplicate names
  (ExperimentSubset #57); column lookups allow for that (A6-16).
- decontPro on a subset uses the subset's main assay only, since a
  subset can't see the main object's altExps (ExperimentSubset #58;
  A6-17).
- `decontXcounts(es, subsetName =)` forwards to
  `assay(x, s, subsetAssayName = "decontXcounts")` (A6-19).
- decontX adds no ExperimentSubset detection code (no `.hasSlot()`, no
  `is(x, "SubsetSingleCellExperiment")` branch) until ExperimentSubset
  #5 adds a shared class. No workaround code for upstream bugs: they are
  filed upstream (ExperimentSubset #2-#5, #57-#59).
- ExperimentSubset support comes last in Phase II, after ExperimentSubset
  #53 is in a release. The Suggests floor is that release; #54 isn't
  needed.

**SCE and decontPro.** decontPro gains `altExpName` and writes its
output into the altExp it read (decided 2026-10-19/20). decontX's SCE
method validates `assayName` against `assayNames(x)` (A6-01, done in
Phase I-b).

**Dependencies.** anndataR and ExperimentSubset go in Suggests, behind
`requireNamespace()` with an error naming the missing package. anndataR
adds nothing beyond decontX's current tree; h5ad input and output need
rhdf5. ExperimentSubset adds 40 packages to the 1.11.2 tree.
SeuratObject stays in Imports for now.

**Tests.** Parity tests per class (same data and seed give results
identical to the matrix method, for `decontX()` and `decontPro()`);
write-back tests for every location in the table; an h5ad round trip
(skipped without rhdf5); Seurat split-layer tests with and without
`batch`; error paths for a missing optional package, an unknown column
name, an unknown `subsetName` and wrong orientation. decontPro fixtures
stay at 10 x 8.

Alternatives considered: `setOldClass("AbstractAnnData")` dispatch
(doesn't dispatch, A6-03); converting through a temporary SCE (extra
copies, loses in-place semantics); ExperimentSubset in Imports (pulls in
TreeSummarizedExperiment and 40 packages); the
`get/setSubsetExperiment()` route (replaced by the direct approach).

## Consequences

- Every supported class runs `decontX()` and `decontPro()` in one call
  and gets its results back in the object.
- decontPro container methods change their return type. Scripts that
  read `decontPro(sce, ...)$decontaminated_counts` break; NEWS marks it
  Breaking, and the vignette has a callout (Phase III-b).
- celda and singleCellTK use the SCE path, which stays bit-identical;
  column-name `z`/`batch` only adds behaviour.
- ExperimentSubset support waits on an upstream release.
- More Suggests packages in CI; tests skip cleanly without them.
- The vignettes' manual round trips are replaced by short "Using your
  object" sections (Phase III-b).

## Open points

- The per-batch UMAP naming for Seurat reductions (sanitized keys) and
  how collisions are handled (A6-11).
- The name of decontPro's output assay or layer for Seurat and AnnData
  (the plan sets `decontPro` for SCE and leaves the final name here).
- Where decontPro's other outputs (ambient and background matrices,
  parameters) live for each class.
- Whether SeuratObject moves to Suggests (A5-09; low priority, about
  110 MB at load).
- The exact ExperimentSubset release used as the Suggests floor.
