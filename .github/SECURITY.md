# Security Policy

## Reporting a vulnerability

Report suspected security issues privately to the maintainer:
Joshua D. Campbell <camp@bu.edu>. Please do not open public GitHub
issues for security problems. You should receive a response within a
week.

## Scope

decontX is an analysis package for single-cell genomics data. The main
risks are around data handling and the build toolchain (compiled C++
and Stan code) rather than network services — it makes no network
connections at runtime.

## Rules for automated tools and AI agents

- Never read, log, or commit credentials, tokens, or private keys.
- Never commit absolute local paths or user-identifying environment
  details.
- Never exfiltrate data files (including test fixtures or user datasets)
  to external services.
- Do not add dependencies or download code/data from unvetted sources;
  dependency changes require an approved ADR.
- Run `/security-review` before releases (see dev/RELEASE.md).
