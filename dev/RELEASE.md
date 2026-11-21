# decontX Release Checklist (Bioconductor)

Bioconductor releases twice a year, ~April and ~October. Check the
current release schedule and **freeze dates** each cycle:
https://bioconductor.org/developers/release-schedule/

## Version scheme

Bioconductor x.y.z: **odd y = devel** (e.g. 1.5.z), **even y = release**
(e.g. 1.4.z). z increments with each change and resets at release. The
y bumps at release time are performed by Bioconductor on their side;
we only bump z in devel.

## Remotes

Development happens on GitHub; Bioconductor's git server is the release
authority. Ensure the `bioc` remote exists:

```
git remote add bioc git@git.bioconductor.org:packages/decontX
```

## Pre-freeze checklist (start ~1 month before the freeze)

1. Sync `devel` with `bioc/devel`; resolve any divergence.
2. Run `make check-full` and `make bioccheck` locally (or inspect the weekly
   CI BiocCheck run).
3. Triage every ERROR/WARNING/NOTE into a fix plan (GitHub issues);
   fix; re-run; land fixes via PRs to `devel`.
4. Run `/security-review` on the release branch diff.
5. Refresh the precomputed full-data articles in `vignettes/articles/`
   (#141, once they exist) on `devel`, so the new release branch carries
   them to the release site (ADR-0007).
6. Update `NEWS.md` with a section for the new release version.
7. Bump `Version:` in DESCRIPTION (z increment) and push `devel` to
   both GitHub and `bioc/devel`.

## At release

- Bioconductor creates the new `RELEASE_X_Y` branch. Add/track it:
  `git fetch bioc && git checkout -b RELEASE_X_Y bioc/RELEASE_X_Y`
- Only bug fixes go to the release branch afterwards (z bumps on both
  branches, cherry-picked or fixed twice — never merge devel into
  release).

## Post-release

- Check the official build report for decontX (release AND devel):
  https://bioconductor.org/checkResults/
  The nightly Bioc build report is the source of truth; CI BiocCheck is
  only the early-warning system.
- Confirm the landing page shows the new version:
  https://bioconductor.org/packages/decontX
- Push the new `RELEASE_X_Y` branch to GitHub if it isn't there yet. The
  push runs the pkgdown workflow, which rebuilds the release site at
  https://campbio.github.io/decontX/ once `bioconductor.org/config.yaml`
  names the new release (re-run it from the Actions tab if the push came
  first). Confirm the root shows the new version and that `/dev/` still
  works.
- Run the dev/AUDIT.md dependency/deprecation audit for the new devel cycle.
