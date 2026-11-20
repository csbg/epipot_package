# Obtain log-fold changes and p-values from a differential epiPot analysis.

Obtain log-fold changes and p-values from a differential epiPot
analysis.

## Usage

``` r
parse_limma_results_contrast(limma_fit, contrast_mat)
```

## Arguments

- limma_fit:

  A fit object returned by limma.

- contrast_mat:

  A contrast matrix.

## Value

A data frame containing, for all cell types (column "cell") and genes
("gene"), log-fold changes ("logFC") and p-values ("P.Value",
"adj.P.Val"), as well as the remaining columns returned by
\`limma::topTable()\`.
