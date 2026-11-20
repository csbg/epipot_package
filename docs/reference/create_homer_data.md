# Save scATAC-seq data for analysis with HOMER.

Save scATAC-seq data for analysis with HOMER.

## Usage

``` r
create_homer_data(mat, limma_res, outpath = "homer_data", min_genes = 200)
```

## Arguments

- mat:

  A matrix with aggregated scATAC-seq counts.

- limma_res:

  Results from the epigenetic potential workflow.

- outpath:

  Path where HOMER input files should be saved.

- min_genes:

  Minimum number of epigenetic potential genes per cell type. HOMER
  input files will not be created if the number is lower.
