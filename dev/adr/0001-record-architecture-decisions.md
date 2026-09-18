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
Architecture Decision Record in `dev/adr/`, using `template.md`
(MADR-style: Context / Decision / Consequences / Status). What requires
one: dependency changes, module/file structure, S4 interface or class
design, Stan model changes, build/deploy machinery (see `README.md` for
the granularity rule). ADRs are proposed via GitHub issues, reviewed
like code, merged via PR, and approved by the maintainer. Numbering is
sequential (0001, 0002, ...), tracked in the README index. Existing ADRs
are immutable history — a change of course gets a new ADR that marks the
old one superseded.

## Consequences

- Structural changes have a required first step (write the ADR), which
  slows them down slightly and deliberately.
- Agents can be pointed at `dev/adr/` for authoritative context; the
  safety rule "no structural refactors without an approved ADR" in
  AGENTS.md is enforceable.
- The directory is part of the repo but excluded from the built package
  via the `^dev$` entry in `.Rbuildignore`.
