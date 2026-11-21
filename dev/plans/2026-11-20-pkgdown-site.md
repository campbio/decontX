# pkgdown site, #84

Design approved in chat, 2026-11-20. The decisions are in ADR-0007
(accepted the same day); this file says what each PR changes and how it
is tested. Three parts:

1. **This branch (`feature/pkgdown-site`), PR "Part of #84":** accept
   ADR-0007, restructure the vignettes, add `_pkgdown.yml` and the URL.
   Version 1.11.38.
2. **`feature/pkgdown-deploy`, PR "Fixes #84":** the deploy workflow and
   the `dev/RELEASE.md` step. Stacked on part 1's branch, because it
   updates the AGENTS.md and roadmap lines part 1 adds; its PR goes into
   `devel` after part 1 merges. No version bump (CI and dev files only).
3. **A release-branch PR into `RELEASE_3_23`:** the same workflow file
   only. No version bump.

Then the maintainer enables Pages and runs the first deploys. The
full-data articles are a separate issue.

No package code changes, so results are unchanged.

## Part 1: vignettes, `_pkgdown.yml`, URL

### Vignettes

- `git mv vignettes/decontX.Rmd vignettes/decontX_RNA.Rmd`. Content
  unchanged except links to the overview where the package is
  introduced; `VignetteIndexEntry` unchanged.
- New `vignettes/decontX.Rmd`, the overview ("Get started"), in
  BiocStyle like the others, following the novice code rules:
  - Introduction: what ambient contamination is (free RNA or ADTs in
    the droplet suspension, captured with each cell) and background
    contamination (non-specific antibody binding, ADTs only).
  - Which method: DecontX for scRNA-seq counts, DecontPro for CITE-seq
    ADT counts; one sentence on each model; what each needs (counts,
    optional cluster or cell-type labels, optional empty droplets).
  - Installation: `BiocManager::install("decontX")`.
  - A quick example: `simulateContamination()` and `decontX()` with
    their defaults (both seeded by their `seed` argument, so no
    `set.seed()`), with the log hidden by the chunk option
    `message = FALSE`, then the estimated contamination plotted against
    the simulated truth. A few lines; runs in seconds.
  - Next steps: links to the two workflow vignettes, by `vignette()`
    call and as relative links that work on Bioconductor and on the
    site.
  - Citing: Yang et al. 2020 (DecontX), Yin et al. 2024 (DecontPro).
  - `sessionInfo()`.
- `vignettes/decontPro.Rmd`: any reference to the RNA vignette points to
  `decontX_RNA`.

### Other references

- `R/decontX-package.R`: `vignette("decontX")` is the overview;
  name `vignette("decontX_RNA")` and `vignette("decontPro")` for the
  workflows. `make docs` regenerates `man/decontX-package.Rd`.
- `README.Rmd` / `README.md`: the vignette list names all three; the
  site link (release and `/dev/`).
- `DESCRIPTION`: `URL: https://campbio.github.io/decontX/,
  https://github.com/campbio/decontX` (site first, as pkgdown expects);
  version 1.11.38.
- `NEWS.md`, under 1.11.38: the RNA workflow vignette is now
  `vignette("decontX_RNA")`; `vignette("decontX")` is a new overview of
  both methods; the package has a documentation site.
- `AGENTS.md`: the layout lists the three vignettes and `_pkgdown.yml`;
  the "no pkgdown site yet" line says the site rules apply and
  `make site-check` runs before PRs that touch the site.
- `dev/ROADMAP.md`: III-c in progress, with the parts above.
- `dev/adr/README.md`: ADR-0007 Accepted, 2026-11-20.

### `_pkgdown.yml`

- `url: https://campbio.github.io/decontX/`
- `template: bootstrap: 5`.
- `navbar`: intro (Get started), reference, articles, news; GitHub link.
- `articles`: two sections, "decontX (RNA)" with `decontX_RNA` and
  "decontPro (Protein)" with `decontPro`. Section titles and
  `navbar:` labels so the Articles menu shows both groups.
- `reference`: the five groups in ADR-0007.
- `.Rbuildignore`: add `^_pkgdown\.yml$` and `^pkgdown$`.

### CI

`lint.yaml` gets a second job, `site-check`, that installs pkgdown and
runs `pkgdown::check_pkgdown()` on pushes and PRs (same triggers as
lint). If `check_pkgdown()` turns out to need the package installed, the
job moves into the Bioconductor docker image instead.

### Tests and checks

No new testthat tests: nothing in `R/` changes. Verification:

- `make site-check` passes (every export indexed, URL consistent).
- `make article FILTER=decontX`, `FILTER=decontX_RNA` and
  `FILTER=decontPro` each render.
- `make docs` changes only `man/decontX-package.Rd` and the README.
- `make test`, `make lint`, `make coverage` (coverage unchanged).
- `make check-full` and `make bioccheck`: no new ERROR, WARNING or NOTE;
  in particular the overview has evaluated code and the vignette count
  is fine.

## Part 2: deploy workflow

`.github/workflows/pkgdown.yaml`, as ADR-0007 describes:

- Triggers: `push` to `devel` and `RELEASE_*`; `workflow_dispatch`.
- `permissions: contents: write` (to push `gh-pages`).
- `concurrency: group: pkgdown, cancel-in-progress: false`.
- A first step picks the mode: `devel` → `devel`; `RELEASE_*` → compare
  with `release_version` from `bioconductor.org/config.yaml` (the same
  download-then-awk as `sync-stable`), `release` if it matches, else
  stop with a notice. Anything else stops.
- The build job runs in `bioconductor/bioconductor_docker:devel` or
  `:<RELEASE_X_Y>`, installs the package and its dependencies through
  BiocManager plus pkgdown, and runs
  `pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)`
  with `PKGDOWN_DEV_MODE` set.
- The build job uploads `docs/` as an artifact. A separate deploy job,
  outside the container (which may lack git or rsync) and the only job
  with `contents: write`, publishes it with
  `JamesIves/github-pages-deploy-action`: development from `docs/dev` to
  `target-folder: dev` (clean inside `dev/` only); release from `docs` to
  the root with `clean-exclude: dev`.
- Actions pinned by major-version tag, as the other workflows are.

`dev/RELEASE.md`:
- Pre-freeze: refresh the precomputed articles on `devel` (once they
  exist).
- Post-release: confirm the push of the new `RELEASE_X_Y` rebuilt the
  root site, and that `/dev/` still works.

Verification before hand-off: `actionlint` if available; the mode logic
checked by running its shell step locally against `devel`,
`RELEASE_3_23`, an old `RELEASE_3_22` and another branch name. The real
deploy can only be tested after the maintainer enables Pages; that is
listed as unverified in the hand-off.

## Part 3: release branch

A `fix/pkgdown-release-deploy` branch from `RELEASE_3_23` with the
workflow file only (and the same file in a later release only if it
changes). PR base `RELEASE_3_23`.
