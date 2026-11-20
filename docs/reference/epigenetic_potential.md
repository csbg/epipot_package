# Calculate the epigenetic potential.

Calculate the epigenetic potential.

## Usage

``` r
epigenetic_potential(
  rna_counts,
  rna_cell_types,
  rna_samples,
  atac_counts,
  atac_cell_types,
  atac_samples,
  verbose = TRUE
)
```

## Arguments

- rna_counts:

  A sparse scRNA-seq count matrix.

- rna_cell_types:

  A vector with cell types whose length is equal to the number of
  columns in \`rna_counts\`.

- rna_samples:

  A vector with sample names whose length is equal to the number of
  columns in \`counts\`.

- atac_counts:

  A sparse scATAC-seq count matrix.

- atac_cell_types:

  A vector with cell types whose length is equal to the number of
  columns in \`atac_counts\`.

- atac_samples:

  A vector with sample names whose length is equal to the number of
  columns in \`atac_counts\`.

- verbose:

  If \`TRUE\`, progress messages are printed.

## Value

A data frame containing, for all cell types (column "cell") and genes
("gene"), log-fold changes ("logFC") and p-values ("P.Value",
"adj.P.Val"), as well as the remaining columns returned by
\`limma::topTable()\`.
