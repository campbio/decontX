<!-- Base branch: devel. Use RELEASE_X_Y only for an approved release fix.
     Never target main/master, which is updated automatically. -->

## What changed and why

## How it was tested

## Checklist

- [ ] Tests added or updated, and `make test` passes
- [ ] `make check-full` and `make bioccheck` pass with no new errors or warnings
- [ ] NEWS.md updated for user-facing changes
- [ ] Version bumped (z) if this will be pushed to Bioconductor
- [ ] Related issue linked
- [ ] ADR linked if this changes structure, dependencies, the S4
      interface, or the Stan model: <!-- dev/adr/NNNN -->
- [ ] **Human judgment**: scientific correctness of any change to the
      DecontX/DecontPro algorithms or their defaults has been reviewed
      by a maintainer (not delegated to tests or agents)
