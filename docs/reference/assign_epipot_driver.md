# Determines the epigenetic potential driver

Determines the epigenetic potential driver

## Usage

``` r
assign_epipot_driver(limma_results, final_results_list)
```

## Arguments

- limma_results:

  contrast results from \`parse_limma_results_contrast()\`.

- final_results_list:

  The RNA list from \`extract_rna_de_overlap_results()\`.

## Value

The same dataframe with an additional column 'Driver'
