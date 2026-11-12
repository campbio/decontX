# decontX docs, part 1 of #83

Approved in chat, 2026-11-11. Audit findings from `dev/audits/2026-10-deep-review.md`
(area 7, and A8-03). Version 1.11.32.

#83 is split in two. This PR (part 1) does the decontX side and the
package-wide docs. Part 2 lands after #81 (decontPro arguments, ADR-0006),
so the decontPro vignette is written once, against the final argument
names and the quiet option: the decontPro vignette rewrite (A7-03 to
A7-05), its FAQ and citation, the rest of A7-11, and A7-23. Part 2 closes
#83; this PR says "Part of #83".

No package code changes, so results are unchanged.

## Changes

### decontX vignette (`vignettes/decontX.Rmd`)

1. **Cell-type mapping in view (A7-01).** Remove the hidden
   `deriveMapping` chunk and the claim that decontX works out the mapping.
   First plot `plotDecontXMarkerPercentage()` by cluster, without
   `groupClusters`, so the reader sees which cluster expresses which
   markers. Then show the code that builds `cellTypeMappings`: one cell
   type at a time, picking the cluster where that type's markers have the
   highest average count. No `lapply()` or closures. (Changed during
   implementation from "detected most often": in PBMC4k, two clusters
   tie at 100% detection of the monocyte markers, which are a large part
   of the ambient RNA, while the average count separates every cell type
   by at least 5-fold.) Cluster numbers can
   change between versions, so the mapping stays computed rather than
   typed in. The paragraph on writing the list by hand stays.
2. **Simpler code (A7-09).** The marker UMAP panel loses its
   `lapply(function(g) ...)` closure: one explicit `plotReducedDim()`
   call per cell type (CD3D, MS4A1, LYZ, GNLY), combined with
   `wrap_plots()`.
3. **One clustering run (A7-06).** The `delta` example reuses the first
   run's clusters, `z = sce$decontX_clusters`, with `embedding = FALSE`.
   That skips the variable genes, PCA and UMAP (the init opt-out from
   #79). The EM still runs a second time, because showing what fixed
   priors do needs one. Both runs then share the same clusters, which
   makes the scatter plot a fairer comparison. The text says why `z` and
   `embedding = FALSE` are there.
4. **"Proportion", not "percentage" (A7-08),** for contamination. The
   percentage of cells that express a marker stays a percentage.
5. **FAQ section (A7-07, A8-03),** before Session information:
   - *Why do my decontaminated counts look higher than the original?*
     (#40) A decontaminated count can't be larger than the observed one,
     because decontX multiplies each count by the probability that it is
     native. After normalization, a cluster's own markers can rise
     because removing contamination shrinks the cell's total. Includes
     the one-line check `max(decontXcounts(sce) - counts(sce))`. The
     decontPro half of #40 (the `plotDensity` bandwidth) goes in part 2.
   - *What does `background` need, and what if my data has one cell
     type?* (#43) Background must be empty droplets only. Genes are
     matched by name: extra genes in the raw matrix are dropped, and
     missing ones give an error. Barcodes that are also in `x` are
     removed from the background. With `background`, a single cluster is
     allowed: `z = rep(1, ncol(sce))`.
   - *Should I run doublet detection before or after decontX?* (#44) The
     text approved on 2026-10-29 (audit, open question 4), unchanged.
   - *How do I run decontX on many samples?* (#45) Pass `batch`: each
     sample is fitted on its own, so the samples don't need to fit in
     one model. Give existing annotations with `z`, and use
     `embedding = FALSE` to skip the UMAP, which is most of the run time
     on large data. Explain how the labels come back: a label used in
     more than one batch is prefixed `"<batch>-<label>"`. If memory is
     still a problem, run the samples one at a time and combine the
     results. I'll claim that this gives the same results as `batch`
     only after a check on simulated data shows it does; otherwise the
     FAQ just suggests it as the fallback.
6. **Citation (A7-12).** A "Citing decontX" section: Yang et al. 2020
   with its DOI, `citation("decontX")`, and one line pointing
   CITE-seq/ADT users to the decontPro vignette and Yin et al. 2024.
7. **Housekeeping (A7-19, #57).**
   - Label every chunk (`session_info`).
   - Wrap prose and code at 80 columns.
   - Fix typos: "what what", "can be plotting", "percetage",
     "decontamined", "suppled", "`groupCluster`" (it's `groupClusters`).
   - The example in the `plotDecontXMarkerPercentage` hints uses
     hard-coded cluster numbers. Reword it so it doesn't name clusters
     that may not match.
   - The Seurat section stays as it is, apart from labels and wrapping,
     since III-b replaces it.

### Package help and roxygen

8. **Package help page** (`R/decontX-package.R`, A7-12): `@references`
   for both papers with DOIs, keeping RStan. One line pointing to the two
   vignettes.
9. **`decontX()` `@return`** (`R/decon.R`, A7-08): `contamination` is
   "the proportion of each cell's counts estimated to be contamination,
   between 0 and 1". The same wording goes in the `plotDecontXContamination()`
   description (`R/plot_decontx.R:3`).
10. **`simulateContamination()`** (A7-10): the `delta` default is
    `c(1, 10)` (the docs say `c(1, 5)`). `@return` lists every element:
    `nativeCounts`, `observedCounts`, `NByC`, `z`, `eta`, `phi`,
    `markers`, `numMarkers` and `contamination`.
11. **`decontXcounts()`** (A7-15, #57): a runnable example. It simulates
    counts, runs `decontX()` with `z = s$z` and `embedding = FALSE` so
    the example is fast, reads the assay with the getter, and sets it
    with the setter.
12. **`retrieveFeatureIndex()`** (A7-18):
    - Fix `\code{link{regex}}` to `\link{regex}`.
    - State the defaults of `by` and `exactMatch`.
    - Make the example deterministic (`set.seed()`, or fixed values).
13. `make docs` regenerates the help pages.

### Package-level files

14. **DESCRIPTION** (A7-21):
    - Title: "Decontamination of Single Cell Genomics Data".
    - Description: add `<doi:10.1186/s13059-020-1950-6>` and
      `<doi:10.1093/nar/gkad1032>` after the two citations.
    - Version 1.11.32.
15. **README** (A7-20): add a short Usage section (simulate counts, run
    `decontX()`, read the results; it must run as pasted) and a Citation
    section with both papers and `citation("decontX")`. Edit
    `README.Rmd`. There's no make target to knit it, and its chunks
    aren't evaluated, so I'll mirror the same text into `README.md` by
    hand. If you'd rather knit `README.md` yourself, say so.
16. **AGENTS.md** (A7-14): DecontPro is Yin et al. 2024, not 2023.
17. **NEWS.md**, under 1.11.32: the vignette shows how to map clusters to
    cell types and gains an FAQ; citations in the vignette, the help page
    and the README; `simulateContamination()` docs give the right `delta`
    default; `decontXcounts()` has an example.
18. **ROADMAP:** record part 1 of #83 under Done; part 2 stays in III-a,
    after #81.

### Maintainer changes at review (2026-11-11)

19. No FAQ entry on decontaminated counts that look higher (#40): one
    report only. This applies to part 2 as well: the decontPro vignette
    gets no #40 entry either.
20. The `background` requirements (empty droplets, barcodes, matching
    genes by name, duplicated names) move from the FAQ into a new
    "With the background matrix" section. "Running decontX" splits into
    "Without the background matrix" and "With the background matrix",
    and the latter shows how to import the raw matrix (the template
    moves there from "Importing data").
21. Multiple samples (#45) becomes its own section under "Other important
    notes", after "Working with Seurat", which now comes first there.
22. The "Note on reproducing results from older versions" is a block
    quote, like the note on matrix input.

## Out of scope

- The decontPro vignette, A7-03 to A7-05, A7-11 and A7-23: part 2, after
  #81.
- The rest of #57: the R dependency and function length.
- The Seurat section rewrite: III-b (#82).

## Tests and checks

This PR changes no package code, so it adds no tests. What it does run:

- `make docs`. Generated files other than `man/` must not change; if
  `src/stanExports_shrinkage.*` changes, restore it.
- `make article FILTER=decontX` to render the vignette. There is no
  `_pkgdown.yml`; if pkgdown can't render without one, I'll ask before
  trying anything else.
- A novice-reader pass on the rendered vignette by a fresh subagent:
  does every chunk run as pasted, and is each step explained?
- `make lint`, `make test`, `make coverage` (unchanged), and before the
  PR, `make check-full` and `make bioccheck`. The new example must run,
  and the decontX vignette must have no chunk-label NOTE.

## Results

decontX, decontPro and simulation results are unchanged. In the
vignette, the `delta` comparison now uses the same clusters in both runs,
so its scatter plot changes slightly. The vignette's other figures
change only where the code changes (the marker UMAP panel and the
cluster barplot).

## Related packages

celda ships an identical copy of `retrieveFeatureIndex()`. The roxygen
fixes need a matching celda issue (drafted for review, not filed).
Nothing else touches celda's forwarding layer.

## Versions

`fix/plot-feature-inputs` (#122, #123) is PR #126, at 1.11.31. This
branch takes 1.11.32, on the assumption that #126 merges first. If it
doesn't, renumber when rebasing.
