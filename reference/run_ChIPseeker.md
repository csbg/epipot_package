# Run the ChIPseeker algorithm

Run the ChIPseeker algorithm

## Usage

``` r
run_ChIPseeker(peaks)
```

## Arguments

- peaks:

  A character vector with peaks given as "chromosome-start-end".

## Value

A data frame with annonated peaks, represented as columns "peak" and
"gene". Peaks from the input for which no annotation could be found are
excluded.
