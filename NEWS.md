# decontX 1.11.1

* The cell cluster initialization used by `decontX()` when `z` is not
supplied (log-normalization, variable gene selection, PCA, and UMAP) has
been migrated from scater/scuttle to the scrapper package, because
scuttle deprecated `logNormCounts()`/`normalizeCounts()` in Bioconductor
3.24. This removes the deprecation warnings from `R CMD check`.
* Results of `decontX()` without `z` may differ slightly from previous
versions (same workflow and parameters; different variable-gene
criterion and UMAP implementation). A new `legacyInit = TRUE` argument
reproduces the previous behavior; it requires the scater package (now in
Suggests) and is kept for backwards compatibility for as long as the
upstream functions remain available. Users supplying their own `z` are
unaffected.
* With the default initialization, the UMAP step itself is deterministic
even when `seed = NULL` (scrapper seeds it explicitly), but overall
`decontX()` results with `seed = NULL` remain non-reproducible because
the EM initialization is unseeded — set a seed for full reproducibility.

# decontX 0.99.5 (2023-10-19)

* First submission to Bioconductor after review. 
