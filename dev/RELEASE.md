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
authority. Ensure the upstream remote exists:

```
git remote add upstream git@git.bioconductor.org:packages/decontX
```

## Pre-freeze checklist (start ~1 month before the freeze)

1. Sync `devel` with `upstream/devel`; resolve any divergence.
2. Run `make check` and `make bioccheck` locally (or inspect the weekly
   CI BiocCheck run).
3. Triage every ERROR/WARNING/NOTE into a fix plan (GitHub issues);
   fix; re-run; land fixes via PRs to `devel`.
4. Run `/security-review` on the release branch diff.
5. Update `NEWS.md` with a section for the new release version.
6. Bump `Version:` in DESCRIPTION (z increment) and push `devel` to
   both GitHub and `upstream/devel`.

## At release

- Bioconductor creates the new `RELEASE_X_Y` branch. Add/track it:
  `git fetch upstream && git checkout -b RELEASE_X_Y upstream/RELEASE_X_Y`
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
- Run the dev/AUDIT.md dependency/deprecation audit for the new devel cycle.
