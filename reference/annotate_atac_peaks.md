# Annotate an ATAC-seq count matrix with gene names.

Annotate an ATAC-seq count matrix with gene names.

## Usage

``` r
annotate_atac_peaks(mat)
```

## Arguments

- mat:

  An ATAC-seq count matrix whose rownames are genomic ranges.

## Value

An ATAC-seq whose rownames are gene names. Rows for which no annotation
could be found are removed.
