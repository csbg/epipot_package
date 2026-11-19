#' Default options for epiPot functions
#'
#' Any function that takes an argument 'options' may override these options.
default_options <- list(
  # conditions
  control_condition = "CTRL",
  disease_conditions = c("AD"),

  # column naming
  separator = ":",

  # contrast patterns (using %s for cell type, %c for control, %d for disease)
  contrast_patterns = list(
    within_control = "%c_ATAC_vs_RNA_%s",
    within_disease = "%d_ATAC_vs_RNA_%s",
    interaction = "%d_vs_%c_%s",
    rna_disease_vs_ctrl = "RNA_%d_vs_%c_%s",
    atac_disease_vs_ctrl = "ATAC_%d_vs_%c_%s"
  ),

  # display names for plots
  display_names = list(
    within_control = "Epigenetic potential - CTRL",
    within_disease = "Epigenetic potential - AD",
    interaction = "Interaction effect",
    rna_disease_vs_ctrl = "Differential expression",
    atac_disease_vs_ctrl = "Differential accessibility",
    rna_included = "Interaction effect DE overlap"
  ),

  # thresholds
  p_value_threshold = 0.05,
  logfc_threshold = 0,

  # wedge definitions
  wedge_logic = list(
    "RNA realized"  = expression(logFC_rna > 0 & abs(logFC_rna) > abs(logFC_atac)),
    "RNA created"   = expression(logFC_rna < 0 & abs(logFC_rna) > abs(logFC_atac)),
    "ATAC realized" = expression(logFC_atac < 0 & abs(logFC_atac) > abs(logFC_rna)),
    "ATAC created"  = expression(logFC_atac > 0 & abs(logFC_atac) > abs(logFC_rna))
  )
)

#' Build contrast names
#'
#' @param pattern_key key of the replacement pattern in the 'contrast_patterns' option
#' @param cell_type name of the cell type
#' @param disease_condition name of the disease condition
#' @param options epiPot options
#'
#' @returns the contrast name
build_contrast_name <- function(pattern_key,
                                cell_type,
                                disease_condition = NULL,
                                options = list()) {
  options <- utils::modifyList(default_options, options)
  pattern <- options$contrast_patterns[[pattern_key]]
  name <- gsub("%s", cell_type, pattern)
  name <- gsub("%c", options$control_condition, name)
  if (!is.null(disease_condition)) {
    name <- gsub("%d", disease_condition, name)
  }
  name
}

#' Create a contrast matrix for a differential epiPot analysis
#'
#' @param model model matrix from `model.matrix()`.
#' @param sample_metadata sample metadata from `create_model_frame()`.
#' @param options epiPot options.
#'
#' @returns a contrast matrix
#' @export
create_contrast_matrix <- function(model, sample_metadata, options = list()) {
  options <- utils::modifyList(default_options, options)

  cells <- unique(sample_metadata$cell)

  make_group <- function(cell, tech, cond) {
    paste(cell, tech, cond, sep = options$separator)
  }

  contrast_matrix <- matrix(
    0,
    nrow = ncol(model),
    ncol = 0,
    dimnames = list(colnames(model), NULL)
  )

  add_contrast <- function(cm, name, terms) {
    v <- stats::setNames(rep(0, nrow(cm)), rownames(cm))
    v[names(terms)] <- terms
    cm <- cbind(cm, v)
    colnames(cm)[ncol(cm)] <- name
    cm
  }

  for (cell in cells) {

    # Control ATAC vs RNA
    g_a_ctrl <- make_group(cell, "ATAC", options$control_condition)
    g_r_ctrl <- make_group(cell, "RNA", options$control_condition)

    if (all(c(g_a_ctrl, g_r_ctrl) %in% colnames(model))) {
      contrast_matrix <- add_contrast(
        contrast_matrix,
        build_contrast_name("within_control", cell, options = options),
        c(stats::setNames(1, g_a_ctrl), stats::setNames(-1, g_r_ctrl))
      )
    }

    # Loop over disease conditions
    for (cond in options$disease_conditions) {

      # Disease ATAC vs RNA
      g_a <- make_group(cell, "ATAC", cond)
      g_r <- make_group(cell, "RNA", cond)

      if (all(c(g_a, g_r) %in% colnames(model))) {
        contrast_matrix <- add_contrast(
          contrast_matrix,
          build_contrast_name("within_disease", cell, cond, options = options),
          c(stats::setNames(1, g_a), stats::setNames(-1, g_r))
        )
      }

      # Disease vs control
      g_a_d <- make_group(cell, "ATAC", cond)
      g_r_d <- make_group(cell, "RNA", cond)

      if (all(c(g_a_d, g_r_d, g_a_ctrl, g_r_ctrl) %in% colnames(model))) {

        # Interaction contrast
        contrast_matrix <- add_contrast(
          contrast_matrix,
          build_contrast_name("interaction", cell, cond, options = options),
          c(stats::setNames(1, g_a_d), stats::setNames(-1, g_r_d),
            stats::setNames(-1, g_a_ctrl), stats::setNames(1, g_r_ctrl))
        )

        # RNA contrast
        contrast_matrix <- add_contrast(
          contrast_matrix,
          build_contrast_name("rna_disease_vs_ctrl", cell, cond, options = options),
          c(stats::setNames(1, g_r_d), stats::setNames(-1, g_r_ctrl))
        )

        # ATAC contrast
        contrast_matrix <- add_contrast(
          contrast_matrix,
          build_contrast_name("atac_disease_vs_ctrl", cell, cond),
          c(stats::setNames(1, g_a_d), stats::setNames(-1, g_a_ctrl))
        )
      }
    }
  }

  contrast_matrix
}

#' Obtain log-fold changes and p-values from a differential epiPot analysis.
#'
#' @param limma_fit A fit object returned by limma.
#' @param contrast_mat A contrast matrix.
#'
#' @return A data frame containing, for all cell types (column "cell") and genes
#'   ("gene"), log-fold changes ("logFC") and p-values ("P.Value", "adj.P.Val"),
#'   as well as the remaining columns returned by `limma::topTable()`.
#' @export
parse_limma_results_contrast <- function(limma_fit, contrast_mat) {
  purrr::map(
    rlang::set_names(colnames(contrast_mat)),
      \(c)
        limma::topTable(limma_fit, coef = c, number = Inf) |>
        tibble::as_tibble()
    ) |>
    purrr::list_rbind(names_to = "coef")
}

#' Extract interaction results
#'
#' @param limma_results Limma results from `parse_limma_results_contrast()`.
#' @param sample_metadata sample metadata from `create_model_frame()`.
#' @param options epiPot options.
#'
#' @returns interaction results
#' @export
#' @importFrom rlang .data
extract_interaction_results <- function(limma_results,
                                        sample_metadata,
                                        options = list()) {
  options <- utils::modifyList(default_options, options)
  final_results_list <- list()

  for (ctx in unique(sample_metadata$cell)) {
    for (cond in options$disease_conditions) {

      interaction_contrast <- build_contrast_name(
        "interaction",
        ctx,
        cond,
        options = options
      )
      control_atac_vs_rna <- build_contrast_name(
        "within_control",
        ctx,
        options = options
      )
      disease_atac_vs_rna <- build_contrast_name(
        "within_disease",
        ctx,
        cond,
        options = options
      )

      required <- c(
        interaction_contrast,
        control_atac_vs_rna,
        disease_atac_vs_rna
      )

      if (!all(required %in% limma_results$coef)) {
        warning("Missing contrasts for ", ctx, " ", cond)
        next
      }

      gg <- dplyr::intersect(
        limma_results |>
          dplyr::filter(
            .data$coef == interaction_contrast,
            .data$adj.P.Val < options$p_value_threshold
          ) |>
          dplyr::pull(.data$genes),
        union(
          limma_results |>
            dplyr::filter(
              .data$coef == control_atac_vs_rna,
              .data$adj.P.Val < options$p_value_threshold,
              .data$logFC > options$logfc_threshold
            ) |>
            dplyr::pull(.data$genes),
          limma_results |>
            dplyr::filter(
              .data$coef == disease_atac_vs_rna,
              .data$adj.P.Val < options$p_value_threshold,
              .data$logFC > options$logfc_threshold
            ) |>
            dplyr::pull(.data$genes)
        )
      )

      final_results <-
        limma_results |>
        dplyr::filter(
          .data$coef == interaction_contrast,
          .data$genes %in% gg
        ) |>
        dplyr::arrange(.data$logFC)

      result_name <- paste(ctx, cond, sep = "|")
      final_results_list[[result_name]] <- final_results
    }
  }

  final_results_list
}

#' Extract interaction results
#'
#' @param limma_results Limma results from `parse_limma_results_contrast()`.
#' @param sample_metadata sample metadata from `create_model_frame()`.
#' @param options epiPot options.
#'
#' @returns interaction results
#' @export
#' @importFrom rlang .data
extract_rna_de_overlap_results <- function(limma_results,
                                           sample_metadata,
                                           options = list()) {
  options <- utils::modifyList(default_options, options)
  final_results_RNA_list <- list()

  for (ctx in unique(sample_metadata$cell)) {
    for (cond in options$disease_conditions) {

      between_rna_contrast <- build_contrast_name(
        "rna_disease_vs_ctrl",
        ctx,
        cond,
        options = options
      )
      interaction_contrast <- build_contrast_name(
        "interaction",
        ctx,
        cond,
        options = options
      )

      if (!all(c(between_rna_contrast, interaction_contrast) %in% limma_results$coef)) {
        next
      }

      gg_RNA <- dplyr::intersect(
        limma_results |>
          dplyr::filter(
            .data$coef == between_rna_contrast,
            .data$adj.P.Val < options$p_value_threshold
          ) |>
          dplyr::pull(.data$genes),
        limma_results |>
          dplyr::filter(
            .data$coef == interaction_contrast,
            .data$adj.P.Val < options$p_value_threshold
          ) |>
          dplyr::pull(.data$genes)
      )

      final_results_RNA <-
        limma_results |>
        dplyr::filter(
          .data$coef == between_rna_contrast,
          .data$genes %in% gg_RNA
        ) |>
        dplyr::arrange(.data$logFC)

      result_name <- paste(ctx, cond, sep = "|")
      final_results_RNA_list[[result_name]] <- final_results_RNA
    }
  }

  final_results_RNA_list
}

#' Make a summary table.
#'
#' @param limma_results Limma results from `parse_limma_results_contrast()`.
#' @param final_results_RNA_list The RNA list from `extract_rna_de_overlap_results()`.
#' @param options epiPot options.
#'
#' @returns A summar table.
#' @export
#' @importFrom rlang .data
create_summary_table <- function(limma_results,
                                 final_results_RNA_list,
                                 options = list()) {
  options <- utils::modifyList(default_options, options)

  purrr::map2(
    options$disease_conditions,
    names(final_results_RNA_list) |>
      stringr::str_split("\\|") |>
      purrr::map_chr(1),
    \(condition, cell_type) {
      c(
        epi_control =
          limma_results |>
          dplyr::filter(
            .data$coef == build_contrast_name(
              "within_control",
              cell_type,
              options = options
            ),
            .data$adj.P.Val < options$p_value_threshold,
            .data$logFC > 0
          ) |>
          nrow(),

        epi_disease =
          limma_results |>
          dplyr::filter(
            .data$coef == build_contrast_name(
              "within_disease",
              cell_type,
              condition,
              options = options
            ),
            .data$adj.P.Val < options$p_value_threshold,
            .data$logFC > 0
          ) |>
          nrow(),

        DE =
          limma_results |>
          dplyr::filter(
            .data$coef == build_contrast_name(
              "rna_disease_vs_ctrl",
              cell_type,
              condition,
              options = options
            ),
            .data$adj.P.Val < options$p_value_threshold,
          ) |>
          nrow(),

        DA =
          limma_results |>
          dplyr::filter(
            .data$coef == build_contrast_name(
              "atac_disease_vs_ctrl",
              cell_type,
              condition,
              options = options
            ),
            .data$adj.P.Val < options$p_value_threshold,
          ) |>
          nrow(),

        interaction =
          limma_results |>
          dplyr::filter(
            .data$coef == build_contrast_name(
              "interaction",
              cell_type,
              condition,
              options = options
            ),
            .data$adj.P.Val < options$p_value_threshold,
          ) |>
          nrow(),

        rna_included =
          final_results_RNA_list[[paste(cell_type, condition, sep = "|")]] |>
          nrow()
      ) |>
        tibble::enframe("contrast", "n_deg") |>
        dplyr::mutate(
          .before = 1,
          condition = condition,
          cell_type = cell_type,
          contrast = dplyr::recode_values(
            .data$contrast,
            "epi_control" ~ options$display_names$within_control,
            "epi_disease" ~ options$display_names$within_disease,
            "DE" ~ options$display_names$rna_disease_vs_ctrl,
            "DA" ~ options$display_names$atac_disease_vs_ctrl,
            "interaction" ~ options$display_names$interaction,
            "rna_included" ~ options$display_names$rna_included
          )
        )
    }) |>
    purrr::list_rbind()
}

#' Determines the epigenetic potential driver
#'
#' @param limma_results contrast results from `parse_limma_results_contrast()`.
#' @param final_results_list The RNA list from `extract_rna_de_overlap_results()`.
#'
#' @returns The same dataframe with an additional column 'Driver'
#' @export
#' @importFrom rlang .data
assign_epipot_driver <- function(limma_results, final_results_list) {
  purrr::imap(
    final_results_list,
    \(tbl, ctx_cond) {
      parts <- stringr::str_split_1(ctx_cond, "\\|")
      cell_type <- parts[1]
      condition <- parts[2]

      df_rna <-
        limma_results |>
        dplyr::filter(
          .data$coef == build_contrast_name(
              "rna_disease_vs_ctrl",
            cell_type,
            condition
          ),
          .data$genes %in% tbl$genes
        ) |>
        dplyr::select(.data$genes, logFC_rna = .data$logFC)

      df_atac <-
        limma_results |>
        dplyr::filter(
          .data$coef == build_contrast_name(
            "atac_disease_vs_ctrl",
            cell_type,
            condition
          ),
          .data$genes %in% tbl$genes
        ) |>
        dplyr::select(.data$genes, logFC_atac = .data$logFC)

      dplyr::inner_join(df_rna, df_atac, by = "genes") |>
        dplyr::mutate(
          .before = 1,
          cell_type = cell_type,
          condition = condition,
          driver = dplyr::case_when(
            logFC_rna > 0 & abs(logFC_rna) > abs(logFC_atac) ~ "RNA realized",
            logFC_rna < 0 & abs(logFC_rna) > abs(logFC_atac) ~ "RNA created",
            logFC_atac < 0 & abs(logFC_atac) > abs(logFC_rna) ~ "ATAC realized",
            logFC_atac > 0 & abs(logFC_atac) > abs(logFC_rna) ~ "ATAC created"
          )
        )
    }
  ) |>
    purrr::list_rbind()
}

#' ggplot2 theme for plotting functions
#'
#' @param base_size font size
#' @param base_family font familty
#' @param grid whether to draw a grid
#' @param rotate whether to rotate x-axis labels
#' @param ... other options passed to the package
#'
#' @returns a ggplot2 theme
theme_epipot <- function(base_size = 10,
                         base_family = "Helvetica",
                         grid = TRUE,
                         rotate = FALSE,
                         ...) {
  new_text_format <- ggplot2::element_text(
    colour = "black",
    size = base_size,
    family = base_family
  )
  res <-
    ggplot2::theme_bw(...) +
    ggplot2::theme(
      axis.text = new_text_format,
      legend.text = new_text_format,
      legend.title = new_text_format,
      strip.text = new_text_format,
      axis.title = new_text_format,
      plot.title = new_text_format
    )

  if (!grid)
    res <- res + ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (rotate)
    res <- res + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

  res
}

#' Plot a summary table
#'
#' @param tbl_summary summary table from `create_summary_table()`
#' @param fill_color bar color
#'
#' @returns a ggplot
#' @export
#' @importFrom rlang .data
plot_summary_table <- function(tbl_summary,
                               fill_color = "#a6cee3") {
  p <-
    ggplot2::ggplot(
      tbl_summary,
      ggplot2::aes(.data$cell_type, .data$n_deg)
    ) +
    ggplot2::geom_col(fill = .data$fill_color) +
    ggplot2::labs(x = NULL, y = "Number of Genes") +
    theme_epipot(
      base_size = 6,
      grid = FALSE,
      rotate = TRUE
    )

  if (dplyr::n_distinct(tbl_summary$condition) > 1) {
    p <-
      p +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$condition, .data$contrast),
        ncol = 3,
      )
  } else {
    p <-
      p +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$contrast),
        ncol = 3
      )
  }

  p
}

#' Plot RNA vs ATAC change for each gene.
#'
#' @param plot_data plot data
#'
#' @returns a ggplot
#' @export
#' @importFrom rlang .data
plot_hexbin_scatter <- function(plot_data) {
  ggplot2::ggplot(plot_data, ggplot2::aes(.data$logFC_rna, .data$logFC_atac)) +
    ggplot2::geom_hex() +
    ggplot2::facet_grid(
      ggplot2::vars(.data$condition),
      ggplot2::vars(.data$cell_type)
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "#e34a33",
      linetype = "solid",
      linewidth = 0.5
    ) +
    ggplot2::geom_abline(
      slope = -1,
      intercept = 0,
      color = "#e34a33",
      linetype = "solid",
      linewidth = 0.5
    ) +
    ggplot2::labs(x = "RNA logFC", y = "ATAC logFC") +
    ggplot2::scale_fill_gradient(
      "Count",
      low = "#bdc9e1",
      high = "#045a8d"
    ) +
    theme_epipot(
      base_size = 6,
      grid = FALSE,
      rotate = TRUE,
    )
}

#' Plot percentages of epigenetic potential drivers
#'
#' @param tbl plot data
#'
#' @returns a ggplot
#' @export
#' @importFrom rlang .data
plot_driver_percentages <- function(tbl) {
  tbl |>
    dplyr::count(.data$cell_type, .data$condition, .data$driver) |>
    dplyr::mutate(
      .by = c(.data$cell_type, .data$condition),
      pct = .data$n / sum(.data$n) * 100
    ) |>
    ggplot2::ggplot(ggplot2::aes(.data$cell_type, .data$pct, fill = .data$driver)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::labs(
      x = NULL,
      y = "Percentage of Genes",
    ) +
    ggplot2::scale_fill_manual(
      "Epigenetic Potential Driver",
      values = c(
        "ATAC created" = "#a6cee3",
        "ATAC realized" = "#b2df8a",
        "RNA created" = "#fb9a99",
        "RNA realized" = "#fdbf6f"
      )
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$condition), ncol = 2) +
    theme_epipot(
      base_size = 6,
      grid = FALSE,
      rotate = TRUE,
    )
}
