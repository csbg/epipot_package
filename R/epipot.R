#' Aggregate a count matrix.
#'
#' @param counts A sparse count matrix.
#' @param cell_types A vector with cell types whose length is equal to the
#'   number of columns in `counts`.
#' @param samples A vector with sample names whose length is equal to the number
#'   of columns in `counts`.
#'
#' @return a sparse matrix with the same number of rows as `counts` and with one
#'   column for each cell type/sample combination. Note that combinations where
#'   all values are zero are removed.
#' @export
#'
#' @examples
aggregate_count_matrix <- function(counts, cell_types, samples) {
  cell_types <- clean_string(cell_types)
  samples <- clean_string(samples)
  colnames(counts) <- paste(cell_types, samples, sep = ":")
  sapply(
    unique(colnames(counts)),
    \(col) Matrix::rowSums(counts[, colnames(counts) == col, drop = FALSE])
  )
}


#' Combine two count matrices.
#'
#' @param rna_mat An RNA-seq count matrix.
#' @param atac_mat An ATAC-seq count matrix.
#'
#' @return A merged matrix that only contains rows whose names appear in both
#'   original matrices.
#' @export
#'
#' @examples
merge_count_matrices <- function(rna_mat, atac_mat) {
  colnames(rna_mat) <- paste0(colnames(rna_mat), ":RNA")
  colnames(atac_mat) <- paste0(colnames(atac_mat), ":ATAC")

  common_genes <- intersect(
    rownames(rna_mat),
    rownames(atac_mat)
  )

  cbind(
    rna_mat[common_genes, ],
    atac_mat[common_genes, ]
  )
}


#' Create a model frame from the merged count matrix.
#'
#' @param mat A count matrix containing RNA- and ATAC-seq data.
#'
#' @return A data frame that can be passed to `model.matrix()`.
#' @export
#'
#' @examples
create_model_frame <- function(mat) {
  tibble::tibble(col = colnames(mat)) |>
    tidyr::separate_wider_delim(
      col,
      delim = ":",
      names = c("cell", "sample", "tech"),
      cols_remove = FALSE
    ) |>
    dplyr::mutate(tech = relevel(factor(tech), ref = "RNA")) |>
    tibble::column_to_rownames("col")
}


#' Obtain log-fold changes and p-values for all genes in all cell types.
#'
#' @param limma_fit A fit object returned by limma.
#'
#' @return A data frame containing, for all cell types (column "cell") and genes
#'   ("gene"), log-fold changes ("logFC") and p-values ("P.Value", "adj.P.Val"),
#'   as well as the remaining columns returned by `limma::topTable()`.
#' @export
#'
#' @examples
parse_limma_results <- function(limma_fit) {
  coefss <- grep(":techATAC", colnames(coef(limma_fit)), value = TRUE)
  purrr::map(
    rlang::set_names(coefss),
    \(c)
      limma::topTable(limma_fit, coef = c, number = Inf) |>
      tibble::as_tibble(rownames = "gene")
  ) |>
    purrr::list_rbind(names_to = "cell") |>
    tidyr::separate_wider_regex(cell, c("cell", cell = ".+", "\\:tech.+"))
}


#' Calculate the epigenetic potential.
#'
#' @param rna_counts A sparse scRNA-seq count matrix.
#' @param rna_cell_types A vector with cell types whose length is equal to the
#'   number of columns in `rna_counts`.
#' @param rna_samples A vector with sample names whose length is equal to the
#'   number of columns in `counts`.
#' @param atac_counts A sparse scATAC-seq count matrix.
#' @param atac_cell_types A vector with cell types whose length is equal to the
#'   number of columns in `atac_counts`.
#' @param atac_samples A vector with sample names whose length is equal to the
#'   number of columns in `atac_counts`.
#' @param verbose If `TRUE`, progress messages are printed.
#'
#' @return A data frame containing, for all cell types (column "cell") and genes
#'   ("gene"), log-fold changes ("logFC") and p-values ("P.Value", "adj.P.Val"),
#'   as well as the remaining columns returned by `limma::topTable()`.
#' @export
#'
#' @examples
epigenetic_potential <- function(rna_counts,
                                 rna_cell_types,
                                 rna_samples,
                                 atac_counts,
                                 atac_cell_types,
                                 atac_samples,
                                 verbose = TRUE) {
  if (verbose)
    msg <- message
  else
    msg <- invisible

  msg("Aggregating count matrices …")
  rna_mat_agg <- aggregate_count_matrix(
    counts = rna_counts,
    cell_types = rna_cell_types,
    samples = rna_samples
  )
  atac_mat_agg <- aggregate_count_matrix(
    counts = atac_counts,
    cell_types = atac_cell_types,
    samples = atac_samples
  )

  msg("Scaling count matrices …")
  rna_mat_z <- scale_rows(rna_mat_agg)
  atac_mat_z <- scale_rows(atac_mat_agg)

  msg("Annotating peaks in ATAC-seq data …")
  atac_mat_z_anno <- annotate_atac_peaks(atac_mat_z)
  merged_mat <- merge_count_matrices(rna_mat_z, atac_mat_z_anno)

  msg("Running linear model …")
  sample_metadata <- create_model_frame(merged_mat)
  model <- model.matrix(
    ~ 0 + cell / tech,
    data = sample_metadata
  )
  fit <-
    merged_mat |>
    lmFit(model) |>
    eBayes()
  res <- parse_limma_results(fit)

  msg("… done!")
  res
}
