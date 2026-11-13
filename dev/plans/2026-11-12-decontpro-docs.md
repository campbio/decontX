# decontPro docs, part 2 of #83

Design approved in chat, 2026-11-12. Audit findings from
`dev/audits/2026-10-deep-review.md` (area 7) and the decontPro rewrite
outline in `dev/audits/2026-10-review-drafts.md`. Version 1.11.34.

Part 1 (#127, 1.11.32) did the decontX side and the package-wide docs.
This part rewrites the decontPro vignette against the arguments from #81
(ADR-0006), finishes A7-11 and A7-23, and closes #83.

No package code changes, so results are unchanged.

## Probe behind the tuning section

On 500 labelled cord_blood droplets (`set.seed(42)`):

| Priors (`deltaSd`, `backgroundSd`) | Ambient share | Background share | Density plots |
|---|---|---|---|
| defaults (2e-5, 2e-6) | 8% | 1% | background peaks of CD3, CD8, CD19 stay |
| 10x (2e-4, 2e-5) | 39% | 6% | those peaks removed; CD4's middle (monocyte) mode stays apart from the lowest |
| 100x (2e-3, 2e-4) | 41% | 17% | over-removal: the CD14 monocyte peak shifts down |

The parts summed to the counts within 4e-12 in every fit. Timings came
from a `load_all()` build (no optimization) and don't count; the audit
measured 7.7 s (defaults) and 5.5 s (larger priors) at 500 droplets on an
installed build.

## Changes

### decontPro vignette (`vignettes/decontPro.Rmd`)

1. **Introduction (A7-04).** The model in plain words: each count is split
   into native (protein on the cell), ambient (free ADTs, one profile
   shared by the dataset) and background (other non-specific binding)
   parts. Cite Yin et al. 2024.
2. **Importing data (A7-05).** `CITEseq("cord_blood", modes =
   c("scADT_Counts", "coldata_scRNAseq"), dry.run = FALSE)`, which
   downloads 0.3 MB instead of 23 MB. Chunk options hide the
   harmonizing message and the `ExperimentList` warning. One object name
   per thing (the old vignette reused `counts` for three objects).
3. **Cell types (A7-05).** The bundled `colData(...)$celltype` labels:
   drop the 158 droplets without a label, then `set.seed()` and sample
   500 droplets, keeping the caution about subsampling. Pass the
   character labels as they are (any labels work since #81). Short
   advice on how fine the labels should be, and that a user's own
   clusters (for example Seurat's `Idents()`) can be passed directly.
   Seurat and dplyr leave this vignette.
4. **Running decontPro (A7-03).** The first fit uses the defaults plus
   `verbose = FALSE`, shown in the code, with one sentence on what it
   hides and that Stan's output is shown anyway if the fit fails.
5. **The output (A7-03).** Each of the four elements:
   `decontaminated_counts`, `ambient_counts` and `background_counts` are
   ADT x droplet matrices with the input's row and column names;
   `parameters` gets a short table of its 8 elements (the two priors,
   `r_est`, `p_est`, `delta_est`, `delta_mean_est`, `background_est`,
   `background_mean_est`). A one-line check that the three parts add up
   to the original counts.
6. **Tuning the priors (A7-04).** `plotDensity()` after the default fit;
   a second fit with larger priors (3x, see Changes at review); its
   `plotDensity()`. The text says which way the priors push (larger
   values let more of a spread-out count be called contamination), uses
   the CD4 three-mode pattern as the QC (only the lowest mode should go),
   and warns that going too far removes real signal, such as part of
   the CD14 monocyte peak. Then `plotBoxByCluster()` by cell type on the
   tuned fit. Two fits in total.
7. **With empty droplets (`background`).** An `eval = FALSE` template
   like decontX's, never run (no empty-droplet demo, decision of
   2026-10-26, confirmed 2026-11-12): import the raw matrix, keep its
   antibody rows, and drop the barcodes of the called cells (decontPro,
   unlike decontX, does not remove them itself). The text says only each
   ADT's total across the empty droplets is used, so `background` sets
   the ambient profile and nothing else, and how rows are matched by
   name. The exact import code is checked against singleCellTK's
   documentation when written.
8. **Runtime and reproducibility.** `seed` and `iter` (defaults 12345
   and 50000): the same inputs and seed give the same result. Runtime
   grows with ADTs x droplets; one measured number for this vignette's
   fit on an installed build. Stan's temporary file is deleted after
   the fit (#46, #81).
9. **FAQ (A7-07).** No #40 entry (maintainer, part 1 review). Candidate
   questions, each kept only if sections 3 and 6 don't already answer
   it: large data (time and memory), using your own clusters, how fine
   cell types should be.
10. **Citing DecontPro (A7-12).** Yin et al. 2024 with
    [doi:10.1093/nar/gkad1032](https://doi.org/10.1093/nar/gkad1032), and
    `citation("decontX")`.
11. **Housekeeping (A7-19).** A descriptive `VignetteIndexEntry`, like
    decontX's; every chunk labelled; prose and code wrapped at 80
    columns; ends with `sessionInfo()`.

### decontPro roxygen (`R/decontPro.R`)

12. **Rest of A7-11.**
    - `background`: only each ADT's total across its droplets (row sums)
      is used, pooled into one ambient profile, so it should hold empty
      droplets only, not the called cells.
    - `@return`: the four elements with their dimensions and dimnames,
      and the 8 elements of `parameters`.
13. **A7-23.** The example uses a 10 x 8 count matrix after `set.seed()`,
    with named arguments and no positional priors: the defaults, unless
    they misbehave on the toy data (then the smallest change that works,
    said in the hand-off). It keeps a `background` matrix of the same
    size.
14. `make docs` regenerates `man/decontPro.Rd`.

### Package-level files

15. **dplyr leaves Suggests** (`DESCRIPTION`). No vignette, test,
    example or R file uses it once the decontPro vignette drops it.
    Normally a dependency change needs an ADR; the maintainer chose to
    record it in this spec instead (2026-11-12). Seurat stays (the
    decontX vignette uses it), and so does SingleCellMultiModal.
16. **Version** 1.11.34.
17. **NEWS.md**, under 1.11.34: the decontPro vignette is rewritten (the
    bundled cell types, every output explained, tuning the priors, a
    `background` template, runtime and reproducibility, FAQ, citation);
    the `decontPro()` help documents every output and what `background`
    is used for, and has a small, seeded example; dplyr is no longer
    suggested.
18. **ROADMAP:** III-a part 2 moves to Done.

## Changes at review (2026-11-12 and 2026-11-13)

Rendering on the installed build showed decontPro running far slower
than the audit measured (134 s per 500-droplet fit, 16 s for the
example; filed as #130). The maintainer then asked whether 500 droplets
change the results. A fit on all 7,700 labelled droplets (about 24
minutes per fit) showed that they do, and that the 10x priors
over-correct:

| 500-droplet fit vs full-data defaults, same droplets | log1p correlation | mean gap in share kept (ADT x cell type, types >= 20 droplets) |
|---|---|---|
| defaults | 0.948 | 0.19 |
| 3x (6e-5, 6e-6) | 0.994 | 0.035 |
| 10x (2e-4, 2e-5) | 0.859 | 0.25 |

On the full data, 10x priors remove the CD4 monocyte peak and shrink
the CD14 high peak. Approved changes:

19. The tuned fit uses 3x priors (`deltaSd = 6e-5`,
    `backgroundSd = 6e-6`), and the text says the result is close to the
    defaults on all droplets (correlation 0.99), with 10x as the
    over-correction example (CD4 middle peak, CD14 high peak).
20. One sentence after the box plot: cell types with few droplets in the
    sample (T-cells 2, Precursors 8) are estimated poorly, another reason
    not to subsample.
21. The runtime section gives the measured times on the current build
    (one to two minutes per 500-droplet fit, about 25 minutes for all
    droplets) until #130 is resolved. NEWS makes no claim about the
    example's run time.
22. The full-data run belongs in the III-c article.

## Out of scope

- `delta_est` and `background_est` have no dimnames. Adding them changes
  the output, so it goes in a follow-up issue if the maintainer agrees.
- Container returns and Seurat or SCE input in the vignette: III-b.
- A decontPro run on real empty droplets: the III-c article (or not at
  all, per the 2026-10-26 decision).

## Tests and checks

No package code changes, so no new tests. What runs:

- `make docs`. Only `man/decontPro.Rd` should change; if
  `src/stanExports_shrinkage.*` changes, restore it.
- The example runs in about a second (timed).
- `make article FILTER=decontPro` renders the vignette, with no Stan
  output and no warnings in the HTML.
- A novice-reader pass on the rendered vignette by a fresh subagent:
  does every chunk run as pasted, and is each step explained?
- `make lint`, `make test`, `make coverage` (unchanged), and before the
  PR, `make check-full` and `make bioccheck`. No chunk-label NOTE.

## Results

decontX, decontPro and simulation results are unchanged. The vignette's
figures change: 500 labelled droplets instead of 1000 Seurat-clustered
ones, a density plot for each of the two fits, and the box plot grouped
by the bundled cell types.

## Related packages

None affected: no exported function, argument or result changes, so celda
and singleCellTK see documentation changes only.
