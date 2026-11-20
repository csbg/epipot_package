# Run the PROGENy analysis.

Run the PROGENy analysis.

## Usage

``` r
run_progeny(counts, limma_results)
```

## Arguments

- counts:

  A merged count matrix as returned by \`merge_count_matrices()\`.

- limma_results:

  Results of the epigenetic potential analysis as returned by
  \`parse_limma_results()\`.

## Value

A list of class \`progeny_results\`.
