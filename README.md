
<!-- README.md is generated from README.Rmd. Please edit that file -->

# decontX

<!-- badges: start -->
<!-- badges: end -->

Methods for decontamination of single cell data

## Installation Instructions

You can install the development version of decontX from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("campbio/decontX")
```

## Use DecontPro

``` r
# Load CITE-seq dataset
library(DropletUtils)
sce_filtered <- read10xCounts('filtered_feature_bc_matrix')
adt_filtered <- sce_filtered[rowData(sce_filtered)$Type == 'Antibody Capture',]

# Cleaning of data as see fit.

# Generate cell clusters
library(Seurat)
library(dplyr)
adt_seurat <- CreateSeuratObject(counts(adt_filtered), assay = 'ADT')
adt_seurat <- NormalizeData(adt_seurat, normalization.method = "CLR", margin = 2) %>%
  ScaleData(assay = "ADT") %>%
  RunPCA(assay = "ADT", features = rownames(adt_seurat),
  reduction.name = "pca_adt", npcs = 10) %>%
  FindNeighbors(dims = 1:10, assay = "ADT", reduction = "pca_adt") %>%
  FindClusters(resolution = 0.2)
  
clusters = as.integer(Idents(adt_seurat))

# Run DecontPro
counts <- as.matrix(counts(adt_filtered))
out <- decontPro(counts, clusters)

decontaminated_counts <- out$decontaminated_counts
```
