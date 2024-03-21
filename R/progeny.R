# Run PROGENy using merged counts matrices and limma results
run_progeny <- function(counts_df, limma_results) {
  progeny_list <- list()
  cell_list <- names(table(limma_results$cell))
  counts <-
    counts_df %>%
    dplyr::mutate_if(~ any(is.na(.x)), ~ if_else(is.na(.x),0,.x)) %>%
    as.matrix()
  net <- get_progeny(organism = 'human', top = 500)
  net$mor <- net$weight
  sample_acts <- run_mlm(mat = counts, net = net, .source = 'source',
                         .target = 'target', .mor = 'weight', minsize = 3)

  # Transform to wide matrix
  sample_acts_mat <-
    sample_acts %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    as.matrix()

  # Scale per feature
  sample_acts_mat <- scale(sample_acts_mat)
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  my_breaks <- c(seq(-3, 0, length.out = ceiling(palette_length/2) + 1),
                 seq(0.05, 3, length.out = floor(palette_length/2)))
  # Plot
  heatmap <-
    pheatmap::pheatmap(sample_acts_mat, border_color = NA,
                       color = my_color, show_rownames = FALSE,
                       breaks = my_breaks)
  progeny_list$heatmap <- heatmap
  non_empty <- c()

  for (cell in cell_list) {
    cells <-  cell #'arterial_endothelial' # cell 'alveolar_type_1'  'airway_smooth_muscle'
    message(cells)
    progeny_list[[cells]] <- list()
    deg <- limma_results %>%
      dplyr::filter(cell == cells & logFC > 0) %>%
      dplyr::select(rn, t) %>%
      dplyr::filter(!is.na(t)) %>%
      column_to_rownames(var = "rn") %>%
      as.matrix()
    if (nrow(deg) == 0) {
      message("This data frame is empty")
    } else if (length(intersect(rownames(deg), net$target)) < 20) {
      message("intersection is less than 20")
    } else {
      message("This data frame is not empty")
      non_empty <- c(non_empty, cells)
      progeny_list[[cells]]$deg <- deg
      # Run mlm
      contrast_acts <- run_mlm(mat = deg, net = net, .source = 'source',
                               .target = 'target', .mor = 'weight', minsize = 5)
      progeny_list[[cells]]$contrast_acts <- contrast_acts
      # Plot
      barplot <-
        ggplot(contrast_acts, aes(x = reorder(source, score), y = score)) +
        geom_bar(aes(fill = score), stat = "identity") +
        scale_fill_gradient2(low = "darkblue", high = "indianred",
                             mid = "whitesmoke", midpoint = 0) +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.x =
                element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
              axis.text.y = element_text(size = 10, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        xlab("Pathways")
      progeny_list[[cells]]$barplot <- barplot

      pathways <- names(table(net$source))
      progeny_list[[cells]]$pathways <- list()


      for (ptw in pathways) {
        pathway <- ptw # 'MAPK'
        df <- net %>%
          dplyr::filter(source == pathway) %>%
          arrange(target) %>%
          mutate(ID = target, color = "3") %>%
          column_to_rownames('target')
        inter <- sort(intersect(rownames(deg),rownames(df)))
        df <- df[inter, ]
        df['t_value'] <- deg[inter, ]
        df <- df %>%
          dplyr::filter(p_value < 0.05) %>%
          mutate(color = if_else(weight > 0 & t_value > 0, '1', color)) %>%
          mutate(color = if_else(weight > 0 & t_value < 0, '2', color)) %>%
          mutate(color = if_else(weight < 0 & t_value > 0, '2', color)) %>%
          mutate(color = if_else(weight < 0 & t_value < 0, '1', color))
        progeny_list[[cells]]$pathways[[pathway]]$df <- df
        weight_plot <-
          ggplot(df, aes(x = weight, y = t_value, color = color)) + geom_point() +
          scale_colour_manual(values = c("red","royalblue3","grey")) +
          geom_label_repel(aes(label = ID)) +
          theme_nb() +
          theme(legend.position = "none",
                panel.border = element_rect(color = "black", fill = NA, size = 0.1),
                plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm")) +
          geom_vline(xintercept = 0, linetype = 'dotted') +
          geom_hline(yintercept = 0, linetype = 'dotted') +
          ggtitle(pathway)
        progeny_list[[cells]]$pathways[[pathway]]$plot <- as.ggplot(weight_plot)
      }
    }
  }
  progeny_list$non_empty <- non_empty

  progeny_list
}

# Save PROGENy volcano plots
save_progeny_volcanos <- function(progeny_results, save_path, tissue){
  path <- save_path # '/media/AGFORTELNY/people/rohit/projects/epigenetic_potential/GSE161383/plots'
  non_empty <- progeny_results[['non_empty']]
  for (cell in non_empty) {
    print(cell)
    plotlist1 <-
      list(progeny_results[[cell]]$pathways$WNT$plot,
           progeny_results[[cell]]$pathways$Hypoxia$plot,
           progeny_results[[cell]]$pathways$PI3K$plot,
           progeny_results[[cell]]$pathways$VEGF$plot,
           progeny_results[[cell]]$pathways$Trail$plot,
           progeny_results[[cell]]$pathways$TNFa$plot,
           progeny_results[[cell]]$pathways$TGFb$plot,
           progeny_results[[cell]]$pathways$EGFR$plot)
    plotlist2 <-
      list(progeny_results[[cell]]$pathways$MAPK$plot,
           progeny_results[[cell]]$pathways$NFkB$plot,
           progeny_results[[cell]]$pathways$p53$plot,
           progeny_results[[cell]]$pathways$Androgen$plot,
           progeny_results[[cell]]$pathways$Estrogen$plot,
           progeny_results[[cell]]$pathways$`JAK-STAT`$plot)
    filedir <- file.path(path, tissue, cell)
    if (!dir.exists(filedir)){
      dir.create(filedir, recursive = TRUE)
    }else{
      print("dir exists")
    }
    filename1 <- file.path(filedir, paste0(cell, '_plotlist1.pdf'))
    title = stringr::str_to_title(gsub('_', ' ', cell))
    filename2 <- file.path(filedir, paste0(cell, '_plotlist2.pdf'))
    ggsave(plot = wrap_plots(plotlist1, ncol = 2) +
             plot_annotation(title, caption = 'PROGENy results',
                             theme=theme(plot.title=element_text(hjust=0.5))) ,
           file=filename1, width = 210, height = 297, units = "mm")
    ggsave(plot = wrap_plots(plotlist2, ncol = 2) ,
           file=filename2, width = 210, height = 297, units = "mm")
  }
}
# Filter cells with no DEGs and process
plot_data <- function(progeny_results){
  non_empty <- progeny_results[['non_empty']]
  plotlist <- list()
  for (cell in non_empty){
    print(cell)
    plotlist[[cell]] <-
      progeny_results[[cell]]$contrast_acts |>
      data.frame() |>
      dplyr::mutate(cells = cell) |>
      dplyr::select(source, cells, score, p_value)
  }
  return(plotlist)
}
# Make dotplots using output from plot_data function
make_dotplot <- function(datalist, title){
  merged_list <- do.call(rbind, datalist)
  heat_mat_row <-
    merged_list |>
    #dplyr::mutate(score = ifelse(p_value < 0.05, score, 0)) |>
    dplyr::select(-p_value) |>
    pivot_wider(names_from = source, values_from = score) |>
    data.frame()
  row.names(heat_mat_row) <- heat_mat_row$cells
  heat_mat_row <- heat_mat_row[,-1]
  clust_row <- hclust(dist(heat_mat_row %>% as.matrix()))
  ddgram_row <- as.dendrogram(clust_row)
  ggtree_plot_row <- ggtree::ggtree(ddgram_row)
  ggtree_plot_row

  heat_mat_col <-
    merged_list |>
    dplyr::mutate(score = ifelse(p_value < 0.05, score, 0)) |>
    dplyr::select(-p_value) |>
    pivot_wider(names_from = source, values_from = score) |>
    data.frame()
  row.names(heat_mat_col) <- heat_mat_col$cells
  heat_mat_col <- heat_mat_col[,-1]
  clust_col <- hclust(dist(heat_mat_col %>% as.matrix()))
  ddgram_col <- as.dendrogram(clust_col)
  ggtree_plot_col <- ggtree::ggtree(ddgram_col)
  ggtree_plot_col

  dotplot <-
    merged_list %>%
    #dplyr::mutate(score = ifelse(p_value < 0.05, score, 0)) |>
    #mutate(cells = factor(cells, levels = clust_row$labels[clust_row$order]) #,
    #source = factor(source, levels = clust_col$labels[clust_col$order]) ) %>%
    ggplot(aes(y=stringr::str_to_title(gsub('_', ' ', cells)),
               x = source, color = score, size = -log10(p_value))) +
    geom_point() +
    geom_point(data=merged_list %>% dplyr::filter(p_value < 0.05),
               pch=21, size= 8, stroke = 1.5, colour="firebrick1") +
    theme_nb() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = title, y = '', x = 'Pathway') +
    theme(axis.ticks = element_blank()) +
    scale_color_gradientn(colours = viridis::viridis(20) ,#viridis::turbo(5),# viridis::viridis(20), #RColorBrewer::brewer.pal(11,"RdBu")
                          limits = c(-2,2),
                          oob = scales::squish,
                          name = 'PROGENy score')
  return(dotplot)

}
# Get gene intersection between progeny set and limma res
progeny_intersection <- function(progeny_res, save_path) {
  require(progeny)
  net <- get_progeny(organism = 'human', top = 500)
  cells <- progeny_res[['non_empty']]
  tissue_list <- list()
  for (cell in cells){
    cell_list <- list()
    degs <- rownames(progeny_res[[cell]]$deg)
    pathways <- names(progeny_res[[cell]]$pathways)
    for (pathway in pathways){
      nets <- net[net$source == pathway,]$target
      inter <- intersect(degs, nets)
      cell_list[[pathway]] <- sort(inter)
    }
    tissue_list[[cell]] <- cell_list
  }
  for (name in names(tissue_list)){
    print(name)
    filedir <- file.path(save_path, name)
    if (!dir.exists(filedir)){
      dir.create(filedir, recursive = TRUE)
    }else{
      print(paste(filedir,"exists"))
    }
    db <- plyr::ldply(tissue_list[[name]], rbind) %>%
      remove_rownames %>%
      column_to_rownames(var=".id")
    filename <- file.path(filedir, paste0(name, '_gene_intersection.csv'))
    fwrite(db, filename)
  }
}
theme_rd <- function(grid = TRUE,
                     rotate = FALSE,
                     ...){
  res <-
    theme_bw(...) +
    theme(
      line = element_line(size = BASE_LINE_SIZE),
      axis.text = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      axis.title = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      legend.background = element_blank(),
      legend.text = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      legend.title = element_text(size = BASE_TEXT_SIZE_PT),
      panel.border = element_rect(size = BASE_LINE_SIZE /2),
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      strip.background = element_blank(),
      strip.text = element_text(
        color = "black",
        size = BASE_TEXT_SIZE_PT
      ),
      strip.text.x = element_text(margin = margin(b = 1, unit = "mm")),
      strip.text.y = element_text(margin = margin(l = 1, unit = "mm"))
    )

  if (!grid)
    res <-
      res +
      theme(panel.grid = element_blank())
  if (rotate)
    res <-
      res +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  res
}
