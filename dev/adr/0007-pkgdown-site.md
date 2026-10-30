# ADR 0007: Publish a pkgdown site on GitHub Pages with release and development versions

## Status

Proposed (2026-10-27). A stub for maintainer review: the decisions below
were made in the plan and during the Phase 0 review (2026-10-18 to
2026-10-27); the open points at the end are not yet decided.

## Context

decontX has no documentation site. Its help and its two Bioconductor
vignettes are the only documentation, and the vignettes must stay short
and fast to build on Bioconductor. The plan
(`dev/plans/2026-10-18-deep-review-and-multi-container-design.md`,
Phase III-c) puts a pkgdown site with full-data articles in scope, so
the vignettes can stay short and link to longer worked examples.

The plan first excluded a development site. On 2026-10-27 the maintainer
added one, as part of moving the documentation to GitHub Pages (plan
section 10).

pkgdown's `development: mode: auto` decides between release and
development from the version number: a fourth component (such as
`.9000`) means development. Bioconductor versions have three parts in
both release and devel (for example 1.10.x and 1.11.x), so `mode: auto`
would call every build a release.

## Decision

We will publish a pkgdown site on GitHub Pages with two versions.

**Sites.**
- The release site, at the root, documents the current Bioconductor
  release. It is rebuilt at each release from the current `RELEASE_*`
  branch.
- The development site, at `/dev/`, is built from `devel` (decided
  2026-10-27).

**Mode by branch.** The build sets the pkgdown development mode from the
branch it builds, not from the version number: `devel` builds go to
`/dev/`, `RELEASE_*` builds to the root.

**Deploy.** A deploy workflow from the shared standards publishes to the
`gh-pages` branch. The `docs/` deny rule in the repo is unaffected.

**Content.**
- `_pkgdown.yml` with a reference index.
- Full-data articles in `vignettes/articles/`: the full PBMC data for
  decontX and the full cord_blood CITE-seq data for decontPro. Heavy
  articles may use precomputed outputs.
- Articles are run and checked locally only. `dev/RELEASE.md` gains a
  step to re-run the articles and rebuild the release site at each
  release.
- The development site renders the short vignettes in CI. The full-data
  articles are never rendered in CI.

**Checks.** CI keeps `pkgdown::check_pkgdown()`. `make site-check` runs
before PRs that touch the site.

**Links.** The site URL goes into DESCRIPTION `URL:` and the README.
After the first deploy, the short vignettes link to the matching
articles. The FAQ stays in the vignettes, because it ships with the
package; articles link back to it.

**Maintainer actions.** Enabling GitHub Pages and the first deploy.

## Consequences

- Users get a browsable reference and full-data examples without
  slowing the Bioconductor build.
- The development site shows unreleased changes, so the site must make
  clear which version a page documents.
- CI gains a site build on `devel`, which renders the short vignettes
  and adds build time to those runs.
- Each release has an extra local step: re-run the articles and rebuild
  the release site.
- `AGENTS.md` notes that there is no pkgdown site; that line and the
  standards' site rules apply once the site exists.

## Open points

- How the full-data articles appear on the development site: copied
  from the last release build, built from precomputed outputs, or
  omitted with a link to the release site.
- The deploy triggers: on every push to `devel`, or on a schedule; and
  whether the release build triggers on a tag, a branch, or by hand.
- How the build learns the current `RELEASE_*` branch when a new
  Bioconductor release is cut.
- The site layout and reference grouping in `_pkgdown.yml`.
