# Package index

## Simple epiPot analysis

- [`epigenetic_potential()`](https://csbg.github.io/epipot_package/reference/epigenetic_potential.md)
  : Calculate the epigenetic potential.
- [`aggregate_count_matrix()`](https://csbg.github.io/epipot_package/reference/aggregate_count_matrix.md)
  : Aggregate a count matrix.
- [`scale_rows()`](https://csbg.github.io/epipot_package/reference/scale_rows.md)
  : Scale the rows of a matrix.
- [`annotate_atac_peaks()`](https://csbg.github.io/epipot_package/reference/annotate_atac_peaks.md)
  : Annotate an ATAC-seq count matrix with gene names.
- [`merge_count_matrices()`](https://csbg.github.io/epipot_package/reference/merge_count_matrices.md)
  : Combine two count matrices.
- [`create_model_frame()`](https://csbg.github.io/epipot_package/reference/create_model_frame.md)
  : Create a model frame from the merged count matrix.
- [`parse_limma_results()`](https://csbg.github.io/epipot_package/reference/parse_limma_results.md)
  : Obtain log-fold changes and p-values for all genes in all cell
  types.

## Differential epiPot analysis

- [`assign_epipot_driver()`](https://csbg.github.io/epipot_package/reference/assign_epipot_driver.md)
  : Determines the epigenetic potential driver
- [`create_contrast_matrix()`](https://csbg.github.io/epipot_package/reference/create_contrast_matrix.md)
  : Create a contrast matrix for a differential epiPot analysis
- [`create_summary_table()`](https://csbg.github.io/epipot_package/reference/create_summary_table.md)
  : Make a summary table.
- [`extract_interaction_results()`](https://csbg.github.io/epipot_package/reference/extract_interaction_results.md)
  : Extract interaction results
- [`extract_rna_de_overlap_results()`](https://csbg.github.io/epipot_package/reference/extract_rna_de_overlap_results.md)
  : Extract interaction results
- [`parse_limma_results_contrast()`](https://csbg.github.io/epipot_package/reference/parse_limma_results_contrast.md)
  : Obtain log-fold changes and p-values from a differential epiPot
  analysis.

## Plotting functions

- [`plot_driver_percentages()`](https://csbg.github.io/epipot_package/reference/plot_driver_percentages.md)
  : Plot percentages of epigenetic potential drivers
- [`plot_hexbin_scatter()`](https://csbg.github.io/epipot_package/reference/plot_hexbin_scatter.md)
  : Plot RNA vs ATAC change for each gene.
- [`plot_summary_table()`](https://csbg.github.io/epipot_package/reference/plot_summary_table.md)
  : Plot a summary table
- [`theme_epipot()`](https://csbg.github.io/epipot_package/reference/theme_epipot.md)
  : ggplot2 theme for plotting functions

## Datasets

- [`rna_data`](https://csbg.github.io/epipot_package/reference/rna_data.md)
  : Sample scRNA-seq dataset
- [`atac_data`](https://csbg.github.io/epipot_package/reference/atac_data.md)
  : Sample scATAC-seq dataset

## Options

- [`default_options`](https://csbg.github.io/epipot_package/reference/default_options.md)
  : Default options for epiPot functions

## Internal functions

- [`build_contrast_name()`](https://csbg.github.io/epipot_package/reference/build_contrast_name.md)
  : Build contrast names
- [`clean_string()`](https://csbg.github.io/epipot_package/reference/clean_string.md)
  : Remove special characters from a string.
- [`granges_to_dataframe()`](https://csbg.github.io/epipot_package/reference/granges_to_dataframe.md)
  : Convert an annotated GRanges object to a data frame.
- [`run_ChIPseeker()`](https://csbg.github.io/epipot_package/reference/run_ChIPseeker.md)
  : Run the ChIPseeker algorithm
