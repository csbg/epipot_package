sanity_plot <- function(matrix = merged_matrix, design = design_mat,
                        gene = gene_name, column_names = pattern,
                        cells = cell_type, ncolumn = facet_cols){
  #gene <- gene_name
  x <- design
  x$E <- as.matrix(matrix)[gene,]
  x$gene <- as.matrix(matrix)[gene,]
  x <- x %>% dplyr::filter(str_detect(x, column_names))
  x$cell <- stringr::str_to_title(gsub("_", " ", x$cell))
  gene_plot <-
    ggplot(x, aes(x=tech, y=gene, fill = tech)) +
    geom_boxplot(lwd=0.3 ,color="black", outlier.colour="firebrick",
                 outlier.shape=8, outlier.size=0.1) + theme_rd()  +
    scale_fill_manual(values = c("grey50","coral")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = cells, y = gene, x = '', fill= "") +
    facet_wrap(.~cell, ncol = ncolumn)
  return(gene_plot)
}
