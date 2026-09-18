# ADR 0001: Record architecture decisions

## Status

Accepted (2026-09-18)

## Context

decontX is maintained by a rotating group of lab members and, increasingly,
AI agents. Decisions about structure (dependencies, file organization, the
S4 interface, the Stan model) were previously made in PR threads and lost.
Agents and new contributors need a durable record of what was decided and
why, so they don't re-litigate settled questions or silently reverse them.

## Decision

We record every structural or dependency decision as a numbered
Architecture Decision Record in `docs/adr/`, using `template.md`
(MADR-style: Context / Decision / Consequences / Status). ADRs are
proposed via GitHub issues, reviewed like code, and merged via PR.
Existing ADRs are immutable history — a change of course gets a new ADR
that supersedes the old one.

## Consequences

- Structural changes have a required first step (write the ADR), which
  slows them down slightly and deliberately.
- Agents can be pointed at `docs/adr/` for authoritative context; the
  safety rule "no structural refactors without an approved ADR" in
  AGENTS.md is enforceable.
- The directory is part of the repo but excluded from the built package
  via `.Rbuildignore`.
