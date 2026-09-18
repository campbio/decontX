# Periodic Dependency / Deprecation Audit

Run once per Bioconductor devel cycle (after each release) or when the
weekly CI BiocCheck job starts failing. This is a *read-only* audit:
findings become GitHub issues (and ADRs where structural) — the audit
itself changes no code.

## Prompt for an agent

> Run `BiocCheck::BiocCheck()` and `BiocManager::valid()` on this
> package. Use r-lib lifecycle practices to find deprecated functions or
> S4 methods and report replacements compatible with the current
> Bioconductor release. Write findings to `agent-log.md`. Do not change
> code — findings become GitHub issues (and ADRs where structural).

## decontX-specific areas to watch

- **rstan / StanHeaders / rstantools**: version bumps regularly break
  compilation of `src/stanExports_shrinkage.cc`; check the rstantools
  changelog and re-generate Stan exports if the scaffolding format
  changed.
- **Matrix**: sparse-matrix class/coercion deprecations (the package
  imports `dgCMatrix` handling); Matrix >= 1.5.3 pinned in DESCRIPTION.
- **Seurat**: major-version API changes affect the `decontPro,Seurat`
  method and `GetAssayData()`-style calls.
- **DelayedArray / SummarizedExperiment / SingleCellExperiment**:
  accessor deprecations surface as new NOTEs in Bioc devel first.
- **celda**: `R/celda_functions.R` mirrors celda internals — check for
  drift between the two packages.

## Output

Append a dated section to `agent-log.md` with: tool versions used,
findings (grouped: errors / warnings / deprecations / drift), and the
list of issues opened.
