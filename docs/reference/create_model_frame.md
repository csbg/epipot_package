# Create a model frame from the merged count matrix.

Create a model frame from the merged count matrix.

## Usage

``` r
create_model_frame(mat, include_condition = FALSE)
```

## Arguments

- mat:

  A count matrix containing RNA- and ATAC-seq data.

- include_condition:

  If true, include the condition as covariate.

## Value

A data frame that can be passed to \`model.matrix()\`.
