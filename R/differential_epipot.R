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
  options <- modifyList(default_options, options)
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
  options <- modifyList(default_options, options)

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
    v <- setNames(rep(0, nrow(cm)), rownames(cm))
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
        c(setNames(1, g_a_ctrl), setNames(-1, g_r_ctrl))
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
          c(setNames(1, g_a), setNames(-1, g_r))
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
          c(setNames(1, g_a_d), setNames(-1, g_r_d),
            setNames(-1, g_a_ctrl), setNames(1, g_r_ctrl))
        )

        # RNA contrast
        contrast_matrix <- add_contrast(
          contrast_matrix,
          build_contrast_name("rna_disease_vs_ctrl", cell, cond, options = options),
          c(setNames(1, g_r_d), setNames(-1, g_r_ctrl))
        )

        # ATAC contrast
        contrast_matrix <- add_contrast(
          contrast_matrix,
          build_contrast_name("atac_disease_vs_ctrl", cell, cond),
          c(setNames(1, g_a_d), setNames(-1, g_a_ctrl))
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
        limma::topTable(fit2, coef = c, number = Inf) |>
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
extract_interaction_results <- function(limma_results,
                                        sample_metadata,
                                        options = list()) {
  options <- modifyList(default_options, options)
  final_results_list <- list()

  for (ctx in unique(sample_metadata$cell)) {
    for (cond in options$disease_conditions) {

      interaction_contrast <- build_contrast_name("interaction", ctx, cond, options = options)
      control_atac_vs_rna <- build_contrast_name("within_control", ctx, options = options)
      disease_atac_vs_rna <- build_contrast_name("within_disease", ctx, cond, options = options)

      required <- c(interaction_contrast, control_atac_vs_rna, disease_atac_vs_rna)

      if (!all(required %in% limma_results$coef)) {
        warning("Missing contrasts for ", ctx, " ", cond)
        next
      }

      gg <- dplyr::intersect(
        dplyr::filter(limma_results,
                      coef == interaction_contrast &
                        adj.P.Val < options$p_value_threshold)$genes,
        union(
          dplyr::filter(limma_results,
                        coef == control_atac_vs_rna &
                          adj.P.Val < options$p_value_threshold &
                          logFC > options$logfc_threshold)$genes,
          dplyr::filter(limma_results,
                        coef == disease_atac_vs_rna &
                          adj.P.Val < options$p_value_threshold &
                          logFC > options$logfc_threshold)$genes
        )
      )

      final_results <- limma_results %>%
        dplyr::filter(coef == interaction_contrast & genes %in% gg) %>%
        arrange(logFC)

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
extract_rna_de_overlap_results <- function(limma_results,
                                           sample_metadata,
                                           options = list()) {
  options <- modifyList(default_options, options)
  final_results_RNA_list <- list()

  for (ctx in unique(sample_metadata$cell)) {
    for (cond in options$disease_conditions) {

      between_rna_contrast <- build_contrast_name("rna_disease_vs_ctrl", ctx, cond, options = options)
      interaction_contrast <- build_contrast_name("interaction", ctx, cond, options = options)

      if (!all(c(between_rna_contrast, interaction_contrast) %in% limma_results$coef)) {
        next
      }

      gg_RNA <- dplyr::intersect(
        dplyr::filter(limma_results,
                      coef == between_rna_contrast &
                        adj.P.Val < options$p_value_threshold)$genes,
        dplyr::filter(limma_results,
                      coef == interaction_contrast &
                        adj.P.Val < options$p_value_threshold)$genes
      )

      final_results_RNA <- limma_results %>%
        dplyr::filter(coef == between_rna_contrast & genes %in% gg_RNA) %>%
        arrange(logFC)

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
create_summary_table <- function(limma_results,
                                 final_results_RNA_list,
                                 options = list()) {
  options <- modifyList(default_options, options)

  map2(
    options$disease_conditions,
    names(final_results_RNA_list) |>
      str_split("\\|") |>
      map_chr(1),
    \(condition, cell_type) {
      c(
        epi_control =
          contrast_res |>
          filter(
            coef == build_contrast_name("within_control", cell_type, options = options),
            adj.P.Val < options$p_value_threshold,
            logFC > 0
          ) |>
          nrow(),

        epi_disease =
          contrast_res |>
          filter(
            coef == build_contrast_name("within_disease", cell_type, condition, options = options),
            adj.P.Val < options$p_value_threshold,
            logFC > 0
          ) |>
          nrow(),

        DE =
          contrast_res |>
          filter(
            coef == build_contrast_name("rna_disease_vs_ctrl", cell_type, condition, options = options),
            adj.P.Val < options$p_value_threshold,
          ) |>
          nrow(),

        DA =
          contrast_res |>
          filter(
            coef == build_contrast_name("atac_disease_vs_ctrl", cell_type, condition, options = options),
            adj.P.Val < options$p_value_threshold,
          ) |>
          nrow(),

        interaction =
          contrast_res |>
          filter(
            coef == build_contrast_name("interaction", cell_type, condition, options = options),
            adj.P.Val < options$p_value_threshold,
          ) |>
          nrow(),

        rna_included = nrow(final_results_RNA_list[[paste(cell_type, condition, sep = "|")]])
      ) |>
        enframe("contrast", "n_deg") |>
        mutate(
          .before = 1,
          condition = condition,
          cell_type = cell_type,
          contrast = recode_values(
            contrast,
            "epi_control" ~ options$display_names$within_control,
            "epi_disease" ~ options$display_names$within_disease,
            "DE" ~ options$display_names$rna_disease_vs_ctrl,
            "DA" ~ options$display_names$atac_disease_vs_ctrl,
            "interaction" ~ options$display_names$interaction,
            "rna_included" ~ options$display_names$rna_included
          )
        )
    }) |>
    list_rbind()
}

#' Determines the epigenetic potential driver
#'
#' @param contrast_res contrast results from `parse_limma_results_contrast()`.
#' @param final_results_list The RNA list from `extract_rna_de_overlap_results()`.
#'
#' @returns The same dataframe with an additional column 'Driver'
#' @export
assign_epipot_driver <- function(contrast_res, final_results_list) {
  imap(
    final_results_list,
    \(tbl, ctx_cond) {
      parts <- str_split_1(ctx_cond, "\\|")
      cell_type <- parts[1]
      condition <- parts[2]

      df_rna <-
        contrast_res |>
        dplyr::filter(
          coef == build_contrast_name("rna_disease_vs_ctrl", cell_type, condition),
          genes %in% tbl$genes
        ) |>
        dplyr::select(genes, logFC_rna = logFC)

      df_atac <-
        contrast_res |>
        dplyr::filter(
          coef == build_contrast_name("atac_disease_vs_ctrl", cell_type, condition),
          genes %in% tbl$genes
        ) |>
        dplyr::select(genes, logFC_atac = logFC)

      inner_join(df_rna, df_atac, by = "genes") |>
        mutate(
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
    list_rbind()
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
  new_text_format <- element_text(
    colour = "black",
    size = base_size,
    family = base_family
  )
  res <-
    theme_bw(...) +
    theme(
      axis.text = new_text_format,
      legend.text = new_text_format,
      legend.title = new_text_format,
      strip.text = new_text_format,
      axis.title = new_text_format,
      plot.title = new_text_format
    )

  if (!grid)
    res <- res + theme(panel.grid = element_blank())
  if (rotate)
    res <- res + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  res
}

#' Plot a summary table
#'
#' @param tbl_summary summary table from `create_summary_table()`
#' @param fill_color bar color
#'
#' @returns a ggplot
#' @export
plot_summary_table <- function(tbl_summary,
                               fill_color = "#a6cee3") {
  p <-
    ggplot(
      tbl_summary,
      aes(cell_type, n_deg)
    ) +
    geom_col(fill = fill_color) +
    labs(x = NULL, y = "Number of Genes") +
    theme_epipot(
      base_size = 6,
      grid = FALSE,
      rotate = TRUE,
      base_family = "Helvetica"
    )

  if (n_distinct(tbl_summary$condition) > 1) {
    p <-
      p +
      facet_wrap(
        vars(condition, contrast),
        ncol = 3,
      )
  } else {
    p <-
      p +
      facet_wrap(
        vars(contrast),
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
plot_hexbin_scatter <- function(plot_data) {
  ggplot(plot_data, aes(logFC_rna, logFC_atac)) +
    geom_hex() +
    facet_grid(vars(condition), vars(cell_type)) +
    geom_abline(
      slope = 1, intercept = 0,
      color = "#e34a33", linetype = "solid", linewidth = 0.5
    ) +
    geom_abline(
      slope = -1, intercept = 0,
      color = "#e34a33", linetype = "solid", linewidth = 0.5
    ) +
    labs(x = "RNA logFC", y = "ATAC logFC") +
    scale_fill_gradient(low = "#bdc9e1", high = "#045a8d", name = "Count") +
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
plot_driver_percentages <- function(tbl) {
  tbl |>
    count(cell_type, condition, driver) |>
    mutate(.by = c(cell_type, condition), pct = n / sum(n) * 100) |>
    ggplot(aes(cell_type, pct, fill = driver)) +
    geom_col(position = "stack") +
    labs(
      x = NULL,
      y = "Percentage of Genes",
    ) +
    scale_fill_manual(
      "Epigenetic Potential Driver",
      values = c(
        "ATAC created" = "#a6cee3",
        "ATAC realized" = "#b2df8a",
        "RNA created" = "#fb9a99",
        "RNA realized" = "#fdbf6f"
      )
    ) +
    facet_wrap(vars(condition), ncol = 2) +
    theme_epipot(
      base_size = 6,
      grid = FALSE,
      rotate = TRUE,
    )
}
