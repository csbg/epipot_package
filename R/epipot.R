

#function to merge scRNA-seq and scATAC-seq merged matrices
create_merged_matrix <- function(atac_pseudobulk_annotated,
                                 rna_seq_pseudobulk_seurat) {
  colnames(rna_seq_pseudobulk_seurat) <-
    paste0(colnames(rna_seq_pseudobulk_seurat),".RNA")
  colnames(atac_pseudobulk_annotated) <-
    paste0(colnames(atac_pseudobulk_annotated),".ATAC")

  merged_sc <- merge(
    rna_seq_pseudobulk_seurat,
    atac_pseudobulk_annotated,
    by = 0
  )

  merged_sc %>%
    remove_rownames %>%
    column_to_rownames("Row.names")
}

create_design_matrix <- function(atac_seq_pseudobulk_annotated,
                                 rna_seq_pseudobulk ){
  atac_seq_pseudobulk_annotated <-
    inject_underscores(atac_seq_pseudobulk_annotated)
  rna_seq_pseudobulk <-
    inject_underscores(rna_seq_pseudobulk)

  design <- rbind(
    data.frame(
      x=names(atac_seq_pseudobulk_annotated),
      row.names = paste0(names(atac_seq_pseudobulk_annotated), ".ATAC"),
      cell = sub("_[^_]+$", "", colnames(atac_seq_pseudobulk_annotated)),#names(atac_seq_pseudobulk_annotated),
      tech="ATAC"),
    data.frame(
      x=names(rna_seq_pseudobulk),
      row.names = paste0(names(rna_seq_pseudobulk),".RNA"),
      cell =  sub("_[^_]+$", "", colnames(rna_seq_pseudobulk)),#names(rna_seq_pseudobulk),
      tech="RNA"
    ) )
  design$tech <- relevel(factor(design$tech), ref="RNA")
  return(design)
}

# Function
parse_limma_res <- function(limma_fit) {
  coefss <- grep(":techATAC", colnames(coef(limma_fit)), value = TRUE)
  coefx <- coefss[1]
  results <- data.table()
  for (coefx in coefss) {
    cell <- gsub("cell(.+)\\:tech.+", "\\1", coefx)
    message(coefx)
    results <-
      rbind(
        results,
        data.table(limma::topTable(limma_fit, coef=coefx, number=Inf),
                   keep.rownames=T,  cell = cell)) #cell_type = ctx,
  }

  results
}

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

# Function
epigenetic_potential <- function(
    atac_count_matrix,
    atac_cell_type_vector,
    atac_sample_vector,
    rna_count_matrix,
    rna_cell_type_vector,
    rna_sample_vector
){
  # Create list object
  epipot <- list()
  # Create aggregate matrix from count matrix              (base R)
  rna_aggr_mat <- matrix_aggregate(
    rna_count_matrix,
    as.character(rna_cell_type_vector),
    as.character(rna_sample_vector))
  atac_aggr_mat <- matrix_aggregate(
    atac_count_matrix,
    as.character(atac_cell_type_vector),
    as.character(atac_sample_vector))
  epipot[["rna_pseudobulk"]] <- rna_aggr_mat
  epipot[["atac_aggr_mat"]] <- atac_aggr_mat
  # Z-score scaling of the aggregate matrices               (base R)
  rna_zr <- data.frame(t(scale(t(rna_aggr_mat))) )
  atac_zr <- data.frame(t(scale(t(atac_aggr_mat))) )
  atac_zr <- atac_zr[grep("\\.\\d+$", row.names(atac_zr), invert = TRUE),]
  # annotate peaks to genes for atac-seq                   (Signac, Repitools, ChIPseeker, EnsDb.Hsapiens.v86)
  atac_zr_anno <- annotate_peaks_pseudobulk(atac_zr)
  # Create merged matrix of scRNA-seq and scATAC-seq        (base R)
  merged_mat <- create_merged_matrix(atac_zr_anno, rna_zr)
  # Create metadata for model matrix                        (base R)
  design_meta <-  create_design_matrix(atac_zr_anno, rna_zr)
  epipot[["design"]] <- design_meta
  # create model matrix
  model <- model.matrix(~ 0 + cell / tech, data=design_meta)
  # Arrange columns as per model matrix                     (base R)
  merged_matrix <- merged_mat[, row.names(model)]
  epipot[["merged_matrix"]] <- merged_matrix
  # Create limma model                                      (limma)
  fit <- lmFit(merged_matrix, model)
  fit <- eBayes(fit)
  # Parse limma model results                               (data.table)
  res <- parse_limma_res(fit)
  epipot[["results"]] <- res
  return(epipot)
}
