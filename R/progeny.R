#' Run the PROGENy analysis.
#'
#' @param counts A merged count matrix as returned by `merge_count_matrices()`.
#' @param limma_results Results of the epigenetic potential analysis as returned
#'   by `parse_limma_results()`.
#'
#' @return A list of class `progeny_results`.
#' @export
run_progeny <- function(counts, limma_results) {
  res <- list()

  message("Analyzing overall dataset")
  counts[is.na(counts)] <- 0
  net <- decoupleR::get_progeny(organism = 'human', top = 500)
  sample_acts <- decoupleR::run_mlm(counts, net, .mor = "weight", minsize = 3)
  res$sample_acts <- sample_acts

  for (cell_type in unique(limma_results$cell)) {
    deg <-
      limma_results |>
      dplyr::filter(cell == cell_type, logFC > 0, !is.na(t)) |>
      dplyr::select(gene, t)
    if (nrow(deg) == 0) {
      message("Skipping ", cell_type, " (no differentially expressed genes)")
    } else if (length(intersect(deg$gene, net$target)) < 20) {
      message("Skipping ", cell_type, " (less than 20 genes in PROGENy database)")
    } else {
      message("Analyzing ", cell_type)
      contrast_acts <-
        deg |>
        tibble::column_to_rownames("gene") |>
        as.matrix() |>
        decoupleR::run_mlm(net, .mor = "weight", minsize = 5)
      pathways <-
        net |>
        dplyr::inner_join(deg, by = dplyr::join_by(target == gene)) |>
        dplyr::filter(p_value < 0.05)
      res[[cell_type]] <- list(
        deg = deg,
        contrast_acts = contrast_acts,
        pathways = pathways
      )
    }
  }

  class(res) <- "progeny_results"
  res
}


#' Plot a PROGENy heatmap.
#'
#' @param x A list as returned by `run_progeny()`.
#'
#' @return A `ComplexHeatmap` object.
#' @export
plot_progeny_heatmap <- function(x) {
  x |>
    tidyr::pivot_wider(
      id_cols = 'condition',
      names_from = 'source',
      values_from = 'score'
    ) |>
    tibble::column_to_rownames('condition') |>
    as.matrix() |>
    scale() |>
    ComplexHeatmap::Heatmap(
      col = circlize::colorRamp2(
        c(-3, 0, 3),
        c("darkblue", "white", "red")
      ),
      show_row_names = FALSE,
      heatmap_legend_param = list(
        title = "score",
        at = c(-3, 0, 3)
      )
    )
}

#' Plot PROGENy pathway scores for a single cell type.
#'
#' @param df A list as returned by `run_progeny()`.
#'
#' @return A `ggplot2` object.
#' @export
plot_progeny_cell_barplot <- function(df) {
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      reorder(source, score),
      score,
      fill = score)
    ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradient2(
      low = "darkblue",
      mid = "whitesmoke",
      high = "indianred"
    ) +
    ggplot2::xlab("pathways") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, size = 10, face = "bold"
      ),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}


#' Plot dots for pathways.
#'
#' @param x A list as returned by `run_progeny()`.
#' @param pathway The pathway to plot.
#'
#' @return A `ggplot2` object.
#' @export
plot_progeny_pathway_dots <- function(x, pathway) {
  x |>
    dplyr::filter(source == {{pathway}}) |>
    dplyr::mutate(
      color = dplyr::case_when(
        weight > 0 & t > 0 ~ '1',
        weight < 0 & t < 0 ~ '1',
        .default = '2'
      )
    ) |>
    ggplot2::ggplot(ggplot2::aes(weight, t, color = color)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
    ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
    ggplot2::geom_point() +
    ggrepel::geom_label_repel(ggplot2::aes(label = target)) +
    ggplot2::scale_colour_manual(values = c("red","royalblue3","grey")) +
    ggplot2::ggtitle(pathway) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}


#' Dotplot of PROGENy results.
#'
#' @param x A list as returned by `run_progeny()`.
#'
#' @return A `ggplot2` object.
#' @export
make_dotplot <- function(x) {
  plot_data <-
    progeny_res[2:length(progeny_res)] |>
    purrr::map(\(x) x$contrast_acts) |>
    purrr::list_rbind(names_to = "cell_type")

  ggplot2::ggplot(plot_data, ggplot2::aes(source, cell_type)) +
    ggplot2::geom_point(ggplot2::aes(color = score, size = -log10(p_value))) +
    ggplot2::geom_point(
      data = plot_data |> dplyr::filter(p_value < 0.05),
      shape = 21,
      size = 8,
      stroke = 1.5,
      color = "firebrick1"
    ) +
    ggplot2::xlab("Pathway") +
    ggplot2::ylab("Cell type") +
    ggplot2::scale_color_viridis_c(
      limits = c(-2, 2),
      oob = scales::squish,
      name = 'PROGENy score'
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.ticks = ggplot2::element_blank()
    )
}
