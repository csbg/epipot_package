# Combine two count matrices.

Combine two count matrices.

## Usage

``` r
merge_count_matrices(rna_mat, atac_mat)
```

## Arguments

- rna_mat:

  An RNA-seq count matrix.

- atac_mat:

  An ATAC-seq count matrix.

## Value

A merged matrix that only contains rows whose names appear in both
original matrices.
