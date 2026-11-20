# epiPot

The `epiPot` package provides tools for analyzing the **epi**genetic
**pot**ential of cells in single-cell data.

## Installation

Install the package from GitHub, e.g.

``` r

remotes::install_github("csbg/epipot_package")
pak::pkg_install("csbg/epipot_package")
renv::install("csbg/epipot_package")
```

## Basic usage

``` r

library(epiPot)

res <- epigenetic_potential(
  rna_counts = rna_data@assays$RNA@counts,
  rna_cell_types = rna_data@meta.data$cells_short,
  rna_samples = rna_data@meta.data$orig.ident,
  atac_counts = atac_data@assays$peaks@counts,
  atac_cell_types = atac_data@meta.data$cell_type,
  atac_samples = atac_data@meta.data$orig.ident
)
```

The vignette describes a
[simple](https://csbg.github.io/epipot_package/articles/epiPot.md)
epiPot analysis in detail.
