# Full-data articles, #141

Design agreed in chat, 2026-11-23. The maintainer asked for the build to
go as far as the knit; they knit the articles themselves until the
shared `make knit-article` target lands
(campbio/r-bioc-dev-standards#21).

## What changes from the issue and ADR-0007

The issue and ADR-0007 named "the full PBMC data" for decontX and "the
full cord_blood data" for decontPro. Probes on 2026-11-23 changed both:

- **decontX:** the RNA vignette already runs on all of pbmc4k, so a
  "full" PBMC article adds nothing. What the vignette can't show is the
  background matrix on real data. Clean 10x data (pbmc4k, the 10x human
  brain 3k nuclei) gives similar answers with and without it. One lung
  nuclei sample from a COVID-19 autopsy (Delorey et al. 2021, GSE171668)
  does not: contamination is high in every cluster, so the profile decontX
  builds from the other clusters isn't clean, and the empty droplets
  change the answer: what decontX removes then matches the empty
  droplets. The article is about that. (A comparison with the CellBender
  output on GEO was dropped on 2026-11-25: it isn't needed to show when
  to use the background matrix.)
- **decontPro:** the vignette stays on cord_blood (500-droplet sample).
  The article uses pbmc_10k_protein_v3 from DropletTestFiles instead of
  all of cord_blood, because it has the raw matrix (spongelets, for
  `background`) and isotype controls, which cord_blood lacks. Its
  cell-type labels come from clustering the ADTs, so the maintainer can
  judge how much work labels are. With coarse labels (one T-cell label
  for CD4 and CD8), DecontPro removes the real peak of CD4 or CD8a
  (which one depends on the label names, #148) and shrinks the naive
  CD45RA peak; the article shows that as a pitfall and names the lost
  marker at knit time.

ADR-0007 gets a dated note saying so, rather than a new ADR.

## Files

All in `vignettes/articles/`:

- `decontX_background.Rmd.orig`, `decontPro_pbmc10k.Rmd.orig`: the
  sources, which run the code.
- `decontX_background.Rmd`, `decontPro_pbmc10k.Rmd`: knitted, results
  only; pkgdown renders these. Produced by the maintainer's knit.
- `decontX_background-figures/`, `decontPro_pbmc10k-figures/`: the
  committed PNGs. Each `.Rmd.orig` sets its own relative `fig.path`, so
  a knit run from `vignettes/articles/` (as the shared target does, and
  as each source's header says for a knit by hand) puts them next to
  the `.Rmd`.

pkgdown copies only the files an article references, so the `.Rmd.orig`
sources don't reach the site.

Other files:

- `.Rbuildignore`: `^vignettes/articles$`.
- `_pkgdown.yml`: `articles/decontX_background` in the "decontX (RNA)"
  group, `articles/decontPro_pbmc10k` in "decontPro (Protein)".
- `vignettes/decontX_RNA.Rmd`, background section, and
  `vignettes/decontPro.Rmd`, "With empty droplets": one sentence linking
  the article. `vignettes/decontX.Rmd`: a sentence naming both articles
  after the list of vignettes. These links are absolute URLs to the
  development site (`https://campbio.github.io/decontX/dev/articles/<name>.html`),
  because the vignettes also render on Bioconductor, and the release
  site has the articles only from the next release on. The articles'
  own links back to the vignettes are relative, so each site links to
  its own copy.
- `dev/adr/0007-pkgdown-site.md`: dated note on the datasets.
- `dev/RELEASE.md` step 5: name the knit command.
- `NEWS.md` and DESCRIPTION: 1.11.40.

## Article mechanics

- Each opens with "Computed with decontX X.Y.Z on <date>", from inline R
  evaluated at knit time, and ends with `sessionInfo()`.
- The extra packages stay out of Suggests: only the local knit needs
  them, and each article lists its own for readers. The lung article
  uses DropletUtils, BiocFileCache, SingleCellExperiment, R.utils and
  ggplot2; the PBMC article DropletTestFiles and DropletUtils. Both
  normalize with scrapper, which decontX imports, rather than with
  hand-written code (changed 2026-11-25).
- Seeds are set before every random step (cell calling, clustering), and
  decontX and decontPro use their own `seed` defaults.
- Cluster names come from the markers, not from cluster numbers, so a
  refresh that renumbers clusters still names them correctly: rules on
  the clusters' mean ADT levels (decontPro), and the highest of scaled
  marker-set scores (decontX). A hidden check stops the decontPro knit
  if the main cell types are no longer found.

## decontX_background: outline

Title: "Using empty droplets with DecontX: single-nucleus lung data".

1. Why the background matrix: DecontX without it builds each cluster's
   contamination from the other clusters; when contamination is high in
   every cluster, or comes from cells or debris not in the cell matrix,
   that profile is off. Links to the vignette for the basics and FAQ.
2. Data: sample 02-P005175-S053-R01 (donor D1, lung, nuclei, 10x v3)
   from GSE171668. The raw matrix downloaded from GEO through
   BiocFileCache, read with `DropletUtils::read10xCounts()`.
3. Calling nuclei: `emptyDropsCellRanger()`; the droplets it doesn't
   call are the empty droplets given as `background`. QC keeps nuclei
   with less than 5% mitochondrial counts; the ones that fail are left
   out of both.
4. decontX without `background`, with its own clusters (clean here: one
   cell type each), then with `background`, reusing the first run's
   clusters as `z` so the empty droplets are the only change. Clusters
   are named by marker-set scores.
5. Comparison: contamination per cell type for each run; what each run
   removes against the empty droplets' profile; detection of
   cell-type markers in other cell types before and after.
6. Summary: when to give decontX the empty droplets.

## decontPro_pbmc10k: outline

Title: "DecontPro on 10x PBMC data with empty droplets".

1. Data: pbmc_10k_protein_v3 (filtered and raw) from DropletTestFiles;
   the 17 ADTs, three of them isotype controls.
2. Cell-type labels: take the ADTs out with `splitAltExps()`, normalize
   them with `scrapper::normalizeAdtCounts.se()` (isotypes left out),
   cluster with scrapper's SNN graph, and label clusters by rules on
   their mean ADT levels.
3. DecontPro with the default ambient estimate.
4. DecontPro with `background`: the spongelets in the raw matrix.
5. Comparison: `plotDensity()`, `plotBoxByCluster()`, share removed per
   ADT; the isotype controls as a check (they should be mostly removed).
6. Pitfall: the same fit with one T-cell label removes the real peak of
   CD4 or CD8a, named at knit time, and shrinks the naive CD45RA peak.
   The `mixed` group shows the same problem for CD16 with the empty
   droplets.
7. Links back to the vignette for the priors and the FAQ.

## Testing and checks

- Before handoff: each `.Rmd.orig`'s code runs top to bottom in a
  scratch directory (purled), so the maintainer's knit doesn't fail on a
  code error, and the text is checked against the results and figures.
- After the maintainer's knit: `make article FILTER=<name>` for both
  articles, and `make site-check` (links resolve, no `.Rmd.orig` on the
  site).
- `make article` for the three edited vignettes, `make check-full`,
  `make bioccheck`: nothing new.

No package code changes, so results are unchanged.
