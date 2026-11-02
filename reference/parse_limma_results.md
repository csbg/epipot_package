# Obtain log-fold changes and p-values for all genes in all cell types.

Obtain log-fold changes and p-values for all genes in all cell types.

## Usage

``` r
parse_limma_results(limma_fit)
```

## Arguments

- limma_fit:

  A fit object returned by limma.

## Value

A data frame containing, for all cell types (column "cell") and genes
("gene"), log-fold changes ("logFC") and p-values ("P.Value",
"adj.P.Val"), as well as the remaining columns returned by
\`limma::topTable()\`.
