# Aggregate a count matrix.

Aggregate a count matrix.

## Usage

``` r
aggregate_count_matrix(counts, cell_types, samples)
```

## Arguments

- counts:

  A sparse count matrix.

- cell_types:

  A vector with cell types whose length is equal to the number of
  columns in \`counts\`.

- samples:

  A vector with sample names whose length is equal to the number of
  columns in \`counts\`.

## Value

a sparse matrix with the same number of rows as \`counts\` and with one
column for each cell type/sample combination. Note that combinations
where all values are zero are removed.
