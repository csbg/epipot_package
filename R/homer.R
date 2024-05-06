#' Save scATAC-seq data for analysis with HOMER.
#'
#' @param mat A matrix with aggregated scATAC-seq counts.
#' @param limma_res Results from the epigenetic potential workflow.
#' @param outpath Path where HOMER input files should be saved.
#' @param min_genes Minimum number of epigenetic potential genes per cell type.
#'   HOMER input files will not be created if the number is lower.
#'
#' @export
create_homer_data <- function(mat,
                              limma_res,
                              outpath = "homer_data",
                              min_genes = 200) {
  annotated_peaks <-
    rownames(mat) |>
    run_ChIPseeker()

  limma_res <-
    limma_res |>
    dplyr::filter(adj.P.Val < 0.05)
  limma_res <- split(limma_res, limma_res$cell)

  dir.create(outpath, showWarnings = FALSE)

  # write HOMER input
  purrr::iwalk(
    limma_res,
    function(df, cell_type) {
      if (nrow(df) < min_genes)
        return
      annotated_peaks |>
        dplyr::filter(gene %in% df$gene) |>
        tidyr::separate_wider_delim(
          peak,
          names = c("chr", "start", "end"),
          delim = "-"
        ) |>
        dplyr::mutate(start = as.integer(start) - 1) |>
        dplyr::select(chr:end) |>
        readr::write_tsv(
          file.path(outpath, paste0(cell_type, ".bed")),
          col_names = FALSE
        )
    }
  )

  # write runscript
  readr::write_file(
    'for file in *.bed; do\n  folder=$(basename "$file" .bed)\n  mkdir -p $folder\n  findMotifsGenome.pl $file hg38 $folder -size 200 -mask\ndone',
    file.path(outpath, "run_homer.sh"),
  )
  invisible()
}


#' Read HOMER results.
#'
#' @param path Pathnwith HOMER results.
#'
#' @return A daraframe with HOMER results.
#' @export
read_homer_results <- function(path = "homer_data") {
  homer_files <- list.files(
    path,
    pattern = "knownResults.txt",
    recursive = TRUE,
    full.names = TRUE
  )
  homer_data <-
    readr::read_tsv(
      homer_files,
      col_names = c("motif", "consensus", "p_val", "log_p", "p_adj",
                    "n_target_seqs", "pct_target_seqs",
                    "n_background_seqs", "pct_background_seqs"),
      skip = 1,
      id = "cell_type",
      show_col_types = FALSE
    ) |>
    dplyr::filter(p_val < 0.05) |>
    tidyr::separate_wider_regex(
      cell_type,
      c(".+/", cell_type = ".+", "/knownResults.txt")
    ) |>
    tidyr::separate_wider_delim(
      motif,
      delim = "/",
      names = c("name_family", "source", "database")
    ) |>
    tidyr::separate_wider_regex(
      name_family,
      pattern = c(name = ".+", "\\(", family = ".+?", "\\)"),
      too_few = "align_start"
    ) |>
    # TODO original approach does not work
    # tidyr::separate_wider_delim(
    #   motif_name,
    #   delim = ",",
    #   names = c("name", "subfamily")
    # ) |>
    # dplyr::mutate(name = str_replace_all(name, "\\(", "_")) |>
    # dplyr::mutate(name = str_replace_all(name, "\\)", "")) |>
    # tidyr::separate(name, sep = '_', c('motif_name','family')) |>
    # tidyr::unite(name, c('motif_name','family'), sep = "_", remove = FALSE) |>
    dplyr::mutate(
      pct_target_seqs = readr::parse_number(pct_target_seqs),
      pct_background_seqs = readr::parse_number(pct_background_seqs),
      odds_ratio = pct_target_seqs / pct_background_seqs
    ) |>
    dplyr::arrange(family)

  homer_data
}


#' Create dotplot of HOMER results.
#'
#' @param motifs A dataframe with HOMER results.
#'
#' @return A `ggplot2` plot.
#' @export
plot_motifs <- function(motifs) {
  plot_data <-
    motifs |>
    dplyr::arrange(family, name) |>
    dplyr::mutate(
      name_family =
        paste0(name, " (", family, ")") |>
        forcats::as_factor()
    )

  ggplot2::ggplot(plot_data, ggplot2::aes(cell_type, name_family)) +
    ggplot2::geom_point(
      ggplot2::aes(color = -log10(p_val), size = odds_ratio)
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(plot_data, p_adj < 0.05),
      shape = 21, size = 4, stroke = 1.5, color = "navy"
    ) +
    khroma::scale_color_hawaii() +
    ggplot2::xlab("cell type") +
    ggplot2::ylab("motif") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}
