# ADR 0009: Raise the scrapper floor to 1.5

## Status

Accepted (2026-11-10)

## Context

ADR-0002 moved the default initialization to scrapper and set the floor
at `scrapper (>= 1.2.0)`. scrapper 1.5 split `runUmap()`'s `seed`
argument into `initialize.seed` and `optimize.seed`, so ADR-0002 kept a
formals check that passes whichever the installed version takes, "until
old versions age out".

They have. Bioconductor 3.23 ships scrapper 1.6.3 and Bioconductor devel
(3.24, where this version of decontX goes) ships 1.7.9. Each decontX
release installs the scrapper of its own Bioconductor release, so no
supported user runs decontX devel against scrapper 1.4 or older. The old
branch can't be reached from any installed scrapper, so tests can't
cover it (A3-11 in `dev/audits/2026-10-deep-review.md`, issue #115).

## Decision

We will raise the floor in DESCRIPTION to `scrapper (>= 1.5.0)` and keep
only the call that passes the seed as `initialize.seed` and
`optimize.seed`. This replaces the last consequence of ADR-0002 (the
formals check); the rest of ADR-0002 stands.

Alternative considered: keeping the check with a test that mocks the old
formals. It would test code no supported install runs, and the floor
would still claim that scrapper 1.2 to 1.4 work.

## Consequences

Results are unchanged: with scrapper >= 1.5 the same arguments reach
`runUmap()` as before. Someone who installs decontX devel into an old
library with scrapper < 1.5 is told to update scrapper at install time,
instead of the code working around it. The release branch
(`RELEASE_3_23`) keeps its own floor. Nothing exported changes, so celda
and singleCellTK are unaffected.
