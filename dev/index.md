# decontX

Methods for decontamination of single cell data. This package implements
both DecontX (Yang et al.,
[2020](https://doi.org/10.1186/s13059-020-1950-6)) for single-cell
RNA-seq data and DecontPro (Yin et al.,
[2024](https://doi.org/10.1093/nar/gkad1032)) for single-cell protein
expression data.

![Overview of decontX. Top: a droplet holds a B cell plus ambient RNA
from other cells, so its counts include T cell and monocyte markers;
decontX() or decontPro() removes them and keeps the B cell marker.
Bottom: on a UMAP before decontamination, contaminated cells blur toward
the middle and a T cell marker shows up in B cells and monocytes; after,
the clusters are tighter and the marker stays in T
cells.](reference/figures/decontX-overview.svg)

## Installation Instructions

You can install the package through
[Bioconductor](https://bioconductor.org/packages/decontX) with:

\
`if`` ``(``!`[`require`](https://rdrr.io/r/base/library.html)`(`[`"BiocManager"`](https://bioconductor.github.io/BiocManager/)`, quietly ``=`` ``TRUE``)``)`\
`    `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"BiocManager"``)`\
\
`BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"decontX"``)`

Alternatively, install the development version from the `devel` branch
on [GitHub](https://github.com/campbio/decontX) using the `devtools`
package. Without `ref = "devel"`,
[`install_github()`](https://devtools.r-lib.org/reference/install-deprecated.html)
installs the current release instead:

\
`# install.packages("devtools")`\
`devtools``::`[`install_github`](https://devtools.r-lib.org/reference/install-deprecated.html)`(``"campbio/decontX"``, ref ``=`` ``"devel"``)`

## Usage

DecontX takes a count matrix or a SingleCellExperiment with genes in
rows and cells in columns. A minimal example on simulated data:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`decontX`](https://campbio.github.io/decontX/)`)`\
\
`# Simulate a contaminated count matrix with 300 cells`\
`sim`` ``<-`` `[`simulateContamination`](https://campbio.github.io/decontX/dev/reference/simulateContamination.md)`(``seed ``=`` ``12345``)`\
\
`# Estimate and remove the contamination`\
`result`` ``<-`` `[`decontX`](https://campbio.github.io/decontX/dev/reference/decontX.md)`(``sim``$``observedCounts``)`\
\
`# The proportion of contamination in each cell, between 0 and 1`\
[`head`](https://rdrr.io/r/utils/head.html)`(``result``$``contamination``)`\
\
`# The decontaminated counts`\
`decontaminated`` ``<-`` ``result``$``decontXcounts`

The vignettes show the full workflow, including a SingleCellExperiment,
cluster labels, empty droplets, and plots of the results.

## Documentation

The package website, <https://campbio.github.io/decontX/>, has the
vignettes and the help pages of every function for the current
Bioconductor release. The development version’s site is at
<https://campbio.github.io/decontX/dev/>.

## Vignettes

Vignettes are available on
[Bioconductor](https://bioconductor.org/packages/decontX):

- Getting started: an overview of DecontX and DecontPro, with a quick
  example.
- decontX (RNA): the DecontX workflow for single-cell RNA-seq data.
- decontPro (Protein): the DecontPro workflow for CITE-seq ADT data.

To build the vignettes when installing from GitHub, add the flag
`build_vignettes = TRUE`:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`devtools`](https://devtools.r-lib.org/)`)`\
[`install_github`](https://devtools.r-lib.org/reference/install-deprecated.html)`(``"campbio/decontX"``, ref ``=`` ``"devel"``, build_vignettes ``=`` ``TRUE``)`

Then vignettes can be accessed through:

\
[`vignette`](https://rdrr.io/r/utils/vignette.html)`(``'decontX'``, package ``=`` ``'decontX'``)`\
[`vignette`](https://rdrr.io/r/utils/vignette.html)`(``'decontX_RNA'``, package ``=`` ``'decontX'``)`\
[`vignette`](https://rdrr.io/r/utils/vignette.html)`(``'decontPro'``, package ``=`` ``'decontX'``)`

## Citation

Please cite the method you used:

- DecontX, for single-cell RNA-seq: Yang S, Corbett SE, Koga Y, Wang Z,
  Johnson WE, Yajima M, Campbell JD (2020). Decontamination of ambient
  RNA in single-cell RNA-seq with DecontX. *Genome Biology*, 21(1), 57.
  <https://doi.org/10.1186/s13059-020-1950-6>
- DecontPro, for single-cell protein expression data: Yin Y, Yajima M,
  Campbell JD (2024). Characterization and decontamination of background
  noise in droplet-based single-cell protein expression data with
  DecontPro. *Nucleic Acids Research*, 52(1), e4.
  <https://doi.org/10.1093/nar/gkad1032>

In R, `citation("decontX")` prints both references.
