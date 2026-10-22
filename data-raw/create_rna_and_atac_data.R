library(Seurat)

DATA_PATH <- "/home/wolfgang/Plus/Projects/epipot_data"

# res <- readRDS(file.path(DATA_PATH, "epipot_results.rds"))
#
# genes_up <-
#   res |>
#   filter(logFC > 0, adj.P.Val <= 0.05) |>
#   slice_max(logFC, n = 5, by = cell, na_rm = TRUE, with_ties = FALSE) |>
#   pull(gene)
#
# genes_down <-
#   res |>
#   filter(logFC < 0, adj.P.Val <= 0.05) |>
#   slice_min(logFC, n = 5, by = cell, na_rm = TRUE, with_ties = FALSE) |>
#   pull(gene)
#
# annotated_peaks <-
#   read_csv("data-raw/annotated_peaks.csv")
#
# annotated_peaks |>
#   filter(gene %in% c(genes_up, genes_down)) |>
#   write_csv("selected_features.csv")


# gene with highest/lowest epigenetic potential per cell type
selected_genes <- read.csv("data-raw/selected_features.csv")$gene
selected_ranges <- read.csv("data-raw/selected_features.csv")$peak

rna_data <-
  readRDS(file.path(DATA_PATH, "GSE161382_rna.rds")) |>
  subset(features = selected_genes)
atac_data <-
  readRDS(file.path(DATA_PATH, "GSE161381_atac.rds")) |>
  subset(features = selected_ranges)

usethis::use_data(rna_data, compress = "xz", overwrite = TRUE)
usethis::use_data(atac_data, compress = "xz", overwrite = TRUE)
