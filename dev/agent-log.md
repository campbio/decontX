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

## 2026-09-21 — Comprehensive engineering review

- **Tools:** R 4.5, BiocCheck (Bioc devel container), `BiocManager::valid()`,
  r-lib lifecycle review, lintr.
- **Scope:** full engineering + targeted method checks (not a pure
  dependency audit); run as sequential subagents per dimension.
- **Outcome:** findings triaged into waves in `dev/ROADMAP.md`; structural
  items flagged as ADR proposals.
- **Full write-up:** `dev/audits/2026-09-comprehensive-review.md`
  (companion drafts in `dev/audits/2026-09-review-drafts.md`).
