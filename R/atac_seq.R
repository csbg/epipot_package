#' Run the ChIPseeker algorithm
#'
#' @param peaks A character vector with peaks given as "chromosome-start-end".
#'
#' @return A data frame with annotated peaks.
run_ChIPseeker <- function(peaks) {
  # select database for ChIPseeker
  edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  ensembldb::seqlevelsStyle(edb) <- "UCSC"

  # convert peaks to genomic ranges
  peak_ranges <-
    tibble::tibble(peaks = peaks) |>
    tidyr::separate(
      peaks,
      into = c("chr", "start", "end"),
      sep = "-"
    ) |>
    GenomicRanges::makeGRangesFromDataFrame()

  # run ChIPseeker
  peaks_annotated <- ChIPseeker::annotatePeak(
    peak_ranges,
    tssRegion = c(-3000, 3000),
    TxDb = edb,
    annoDb = "org.Hs.eg.db",
    verbose = FALSE
  )

  # convert ChIPseeker results to data frame
  Repitools::annoGR2DF(peaks_annotated@anno) |>
    tibble::as_tibble() |>
    tidyr::unite(
      chr,
      start,
      end,
      col = "peak",
      sep = "-",
      remove = FALSE
    )
}

#' Annotate an ATAC-seq count matrix with gene names.
#'
#' @param mat An ATAC-seq count matrix whose rownames are genomic ranges.
#'
#' @return An ATAC-seq whose rownames are gene names. Rows for which no
#'   annotation could be found are removed.
#' @export
#'
#' @examples
annotate_atac_peaks <- function(mat) {
  # get ChIPseeker annotations from genomic ranges stored in rownames
  peak_annotations <-
    rownames(mat) |>
    run_ChIPseeker()

  # for each genomic range, reduce annotations to a single gene
  # by selecting the one with the smallest distance to TSS
  annotated_peaks <-
    peak_annotations |>
    dplyr::filter(!is.na(SYMBOL)) |>
    dplyr::arrange(SYMBOL, abs(distanceToTSS)) |>
    dplyr::distinct(SYMBOL, .keep_all = TRUE) |>
    dplyr::select(peak, SYMBOL) |>
    dplyr::arrange(peak)

  # filter count matrix so it only contains rows
  # for which an annotation was found,
  # and replace rownames with these genes
  res <- mat[annotated_peaks$peak, ]
  rownames(res) <- annotated_peaks$SYMBOL
  res
}
