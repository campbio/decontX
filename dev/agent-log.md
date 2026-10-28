# Agent Log

Append-only log of automated / agent-run maintenance passes on decontX —
primarily the periodic dependency & deprecation audit described in
`dev/AUDIT.md`. Newest entry first.

Each entry records: the date, the tool versions used, findings grouped as
**errors / warnings / deprecations / drift**, and the list of issues (and
ADRs) opened. Full write-ups of large reviews live in `dev/audits/`; this
file is the running index of when each pass happened and what came out of
it.

<!-- Template for a new entry (copy above the previous one):

## YYYY-MM-DD — <short title>

- **Tools:** R x.y.z, BiocCheck x.y.z, BiocManager::valid(), lifecycle ...
- **Errors:** ...
- **Warnings:** ...
- **Deprecations:** ...
- **Drift:** ...
- **Issues/ADRs opened:** #NN, ADR-XXXX ...
- **Full write-up:** dev/audits/YYYY-MM-...md (if any)

-->

## 2026-10-27 — Deep review (Phase 0)

- **Tools:** R 4.6.1, local Bioc 3.23 (release) library, BiocCheck
  1.48.1, covr, lintr, bench/profvis via `dev/bench/bench-decontx.R`
  (installed -O2 builds); ExperimentSubset 1.23.16 (PR #54 branch),
  anndataR 1.2.2, SeuratObject 5.4.0 in scratch libraries.
- **Scope:** the Phase 0 review of the 2026-10-18 plan: baseline, status
  of the 2026-09 findings, and eight areas (R bugs, compiled code,
  tests, code quality, speed and memory, object coupling, docs,
  issue-driven methods checks), run as sequential subagents.
- **Errors:** none (`make check` OK, devel CI green).
- **Warnings:** two test warnings (#55, Seurat coercion); BiocCheck
  ERROR and WARNING are environmental.
- **Deprecations:** `Eigen::MappedSparseMatrix` (D2-10); scater
  normalization functions used by `legacyInit` (warnings shown, no
  change).
- **Drift:** local library is Bioc 3.23 while the package targets 3.24;
  a devel library is set up before Phase I.
- **Outcome:** no P0; 11 P1 bugs, mostly silent input-handling failures;
  coverage 74%; initialization is about 90% of large runs.
- **Issues/ADRs opened:** campbio/ExperimentSubset #2-#5 (2026-10-19)
  and #57-#59 (2026-10-27). decontX issues, replies and ADR-0004 to
  ADR-0007 are proposed in the audit, which the maintainer approved on
  2026-10-27 with a pkgdown development site added.
- **Full write-up:** `dev/audits/2026-10-deep-review.md` (evidence in
  `dev/audits/2026-10-review-drafts.md`); `dev/ROADMAP.md` rewritten.

## 2026-09-21 — Comprehensive engineering review

- **Tools:** R 4.5, BiocCheck (Bioc devel container), `BiocManager::valid()`,
  r-lib lifecycle review, lintr.
- **Scope:** full engineering + targeted method checks (not a pure
  dependency audit); run as sequential subagents per dimension.
- **Outcome:** findings triaged into waves in `dev/ROADMAP.md`; structural
  items flagged as ADR proposals.
- **Full write-up:** `dev/audits/2026-09-comprehensive-review.md`
  (companion drafts in `dev/audits/2026-09-review-drafts.md`).
