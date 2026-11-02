# Scale the rows of a matrix.

Scale the rows of a matrix.

## Usage

``` r
scale_rows(m)
```

## Arguments

- m:

  A numeric matrix

## Value

A matrix whose rows have been scaled.

## Examples

``` r
scale_rows(diag(4))
#>      [,1] [,2] [,3] [,4]
#> [1,]  1.5 -0.5 -0.5 -0.5
#> [2,] -0.5  1.5 -0.5 -0.5
#> [3,] -0.5 -0.5  1.5 -0.5
#> [4,] -0.5 -0.5 -0.5  1.5
#> attr(,"scaled:center")
#> [1] 0.25 0.25 0.25 0.25
#> attr(,"scaled:scale")
#> [1] 0.5 0.5 0.5 0.5
```
