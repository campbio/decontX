# ADR 0007: Publish a pkgdown site on GitHub Pages with release and development versions

## Status

Accepted (2026-11-20). Proposed 2026-10-27 as a stub from the deep
review; the maintainer settled its open points on 2026-11-20 (#84).

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

The shared standards (r-bioc-dev-standards) have no deploy workflow, so
decontX carries its own. The current release branch, `RELEASE_3_23`
(1.10.0), has no `_pkgdown.yml` and the old vignette layout.

## Decision

We will publish a pkgdown site on GitHub Pages, at
`https://campbio.github.io/decontX/`, with two versions. The lab's
`https://www.camplab.net/decontx/` will point to it (maintainer action).

### Sites and mode

- The **release site**, at the root, documents the current Bioconductor
  release. It is built from the current `RELEASE_*` branch.
- The **development site**, at `/dev/`, is built from `devel`.
- The build sets the mode from the branch, with the `PKGDOWN_DEV_MODE`
  environment variable: `devel` for `devel` builds, which pkgdown writes
  under `dev/`, and `release` for `RELEASE_*` builds. `_pkgdown.yml` has
  no `development:` mode.

### Deploy workflow

One workflow, `.github/workflows/pkgdown.yaml`, in decontX:

- **Development build:** on every push to `devel`, and on demand
  (`workflow_dispatch`). It deploys only into `dev/` on the `gh-pages`
  branch and never touches the root.
- **Release build:** on every push to a `RELEASE_*` branch, and on
  demand. It first reads `release_version` from
  `https://bioconductor.org/config.yaml`, as the shared `sync-stable`
  workflow does, and stops unless the branch is the current release, so
  a late fix on an old release branch can't overwrite the root. It
  deploys to the root and keeps `dev/`.
- When Bioconductor cuts a new release, pushing the new `RELEASE_X_Y`
  branch to GitHub rebuilds the root; no extra step.
- Each build runs in the Bioconductor docker image for its branch
  (`devel`, or the image tagged with the release branch name), as
  R-CMD-check does.
- Builds are grouped by branch and deploys by site, so a push to one
  branch never cancels the other's waiting build or deploy. The two
  sites' deploys write to separate folders of `gh-pages`, and the deploy
  action retries a rejected push.
- No scheduled builds: nothing changes between pushes.
- Pull requests never deploy. Their CI runs `pkgdown::check_pkgdown()`
  only, and `make site-check` runs before PRs that touch the site.
- The same workflow file is added to `RELEASE_3_23` in its own commit
  (no version bump, as for 7e14846). With no `_pkgdown.yml` there,
  pkgdown builds a default site from the 1.10.0 README, help and
  vignettes, so the root works from the first deploy. The full layout
  reaches the root with the first release cut after this lands (Bioc
  3.24).

Alternatives: building the release site locally and pushing it by hand
(a second deploy path that depends on one machine); a shared deploy
workflow in r-bioc-dev-standards (worth doing once a second package
needs it, not before).

### Content and layout

- **Navbar:** Get started, Reference, Articles, Changelog (from
  `NEWS.md`), and a GitHub link. Bootstrap 5, otherwise pkgdown's
  defaults.
- **Get started** is a short `vignettes/decontX.Rmd`: what ambient and
  background contamination are, when to use DecontX (scRNA-seq) and
  when DecontPro (CITE-seq ADT), installation, links to the two
  workflow vignettes, how to cite each method, and a tiny simulated
  example ending with `sessionInfo()` (so BiocCheck has evaluated code).
- The **RNA workflow vignette** moves to `vignettes/decontX_RNA.Rmd`. It
  stays a package vignette, so it still ships with decontX and appears
  on Bioconductor; `vignette("decontX")` opens the overview instead.
- **Articles menu**, two groups:
  - "decontX (RNA)": `decontX_RNA`, and later the full-data PBMC
    article.
  - "decontPro (Protein)": `decontPro`, and later the full-data
    cord_blood article.
- **Reference index:**

  | Group | Topics |
  |---|---|
  | Decontaminate RNA (DecontX) | `decontX`, `decontXcounts`, `simulateContamination` |
  | Plot DecontX results | `plotDecontXContamination`, `plotDecontXMarkerExpression`, `plotDecontXMarkerPercentage` |
  | Decontaminate protein (DecontPro) | `decontPro`, `plotDensity`, `plotBoxByCluster` |
  | Utilities | `retrieveFeatureIndex` |
  | Package | `decontX-package` |

  `decontSpa()` gets its own "Spatial" group when it lands.

Alternative: moving the vignettes into `vignettes/articles/`. Rejected:
that folder is left out of the package build, so they would no longer
ship with decontX or appear on Bioconductor, `vignette()` would not find
them, BiocCheck needs a real vignette, and the FAQ would leave the
package.

### Full-data articles

- Full-data articles live in `vignettes/articles/` (left out of the
  build): the full PBMC data for decontX and the full cord_blood
  CITE-seq data for decontPro. They are fuller versions of the two
  vignettes; they are their own issue.
- They are **precomputed**. Each source is
  `vignettes/articles/<name>.Rmd.orig`; it is knitted locally into a
  plain `vignettes/articles/<name>.Rmd` that already holds the results,
  with its figures committed beside it. CI renders both sites from these
  files without computing anything, so the development and release
  sites show the same articles. Each article states the decontX version
  that produced it.
- The articles are refreshed on `devel` before each release freeze, so
  the new release branch carries them. `dev/RELEASE.md` gets that step.
  The knit runs through a make target that the maintainer adds to the
  Makefile with the articles.
- *Note, 2026-11-23 (#141):* the datasets changed before the articles
  were written. The RNA vignette already runs on all of pbmc4k, so the
  decontX article is about the `background` matrix instead, on one
  single-nucleus lung sample from GSE171668, where the empty droplets
  change the results. The decontPro article uses pbmc_10k_protein_v3
  from DropletTestFiles instead of all of cord_blood, for its raw
  matrix and isotype controls. The knit target is a shared one,
  `make knit-article FILTER=<name>`, from the standards
  (campbio/r-bioc-dev-standards#21), rather than one in this Makefile.
  Details in `dev/plans/2026-11-23-full-data-articles.md`.
- Alternatives: rendering the articles in CI (too slow, and needs the
  full data downloaded on every build); leaving them off the development
  site (the release build would still need them precomputed or built
  locally, so nothing is saved).

### Links

The site URL goes into DESCRIPTION `URL:`, the README and the `url:`
field of `_pkgdown.yml`, with the vignette restructure. The vignettes
link to the full-data articles once those exist. The FAQ stays in the
vignettes, because it ships with the package; articles link back to it.

### Maintainer actions

Enabling GitHub Pages (from the `gh-pages` branch), the first deploys,
and pointing camplab.net/decontx at the site.

## Consequences

- Users get a browsable reference and full-data examples without
  slowing the Bioconductor build.
- The development site shows unreleased changes; pkgdown's development
  mode marks its version in the navbar.
- CI gains a site build on each push to `devel` (a few minutes: the
  Stan compile and the RNA vignette dominate) and on release-branch
  pushes.
- Until the first deploy, DESCRIPTION's URL points to a page that does
  not exist yet; the vignette restructure and the deploy workflow should
  land close together.
- `vignette("decontX")` opens the overview; old references to the RNA
  workflow by that name need to change.
- Each release has an extra local step: refresh the precomputed
  articles before the freeze. Their figures add a few hundred KB to the
  repository per refresh, but nothing to the package tarball.
- `AGENTS.md` notes that there is no pkgdown site; that line changes, and
  the standards' site rules apply, once `_pkgdown.yml` exists.
