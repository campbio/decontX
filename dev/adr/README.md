# Architecture Decision Records

This directory is the durable record of structural decisions in decontX.
It is append-only: decisions are never edited after acceptance — a change
of course is a new ADR that marks the old one superseded.

## When an ADR is required

ADRs are for decisions that are **hard to reverse** or that **future
contributors will question**: dependency changes (anything touching
DESCRIPTION Imports/LinkingTo), module/file structure, S4 interface or
class design, changes to the Stan model, and build/deploy machinery.
Routine choices (naming a helper, fixing a bug, adding a test) do not
get ADRs.

Process: propose via a GitHub issue, draft from `template.md` (the
`adr-author` skill helps), land via PR. The maintainer approves.
Numbering is sequential, zero-padded to four digits.

## Index

| # | Title | Status | Date |
|---|---|---|---|
| [0001](0001-record-architecture-decisions.md) | Record architecture decisions | Accepted | 2026-09-18 |
| [0002](0002-migrate-initialization-to-scrapper.md) | Migrate cluster initialization from scater/scuttle to scrapper | Accepted | 2026-09-21 |
