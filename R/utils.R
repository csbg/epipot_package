clean_vector <- function(vector) {
  vector <- gsub(" ", "_", vector)
  vector <- gsub("/", "_", vector)
  vector <- gsub("\\-", "_", vector)
  vector <- gsub("\\(", "", vector)
  vector <- gsub("\\)", "", vector)
  vector <- gsub(",", "", vector)
  vector
}

inject_underscores <- function(df) {
  colnames(df) <- gsub("\\.", "_", colnames(df))
  df
}

# Function to aggregate scRNA/scATAC count matrix
# result will be a n_features by (n_cell_types * n_samples) data frame
# TODO: some cell type/sample combinations are removed
matrix_aggregate <- function(counts, cell_types, samples) {
  mat <- counts
  samples <- as.numeric(factor(samples))
  cell_types <- as.character(clean_vector(cell_types))
  samples <- as.character(clean_vector(samples))
  colnames(mat) <- paste(cell_types, samples, sep = "_")
  mat_aggr <-
    sapply(
      unique(colnames(mat)),
      function(x, mat) rowSums(mat[,colnames(mat) == x, drop=FALSE]),
      mat = mat)
  count_aggr_df <- data.frame(mat_aggr)
  return(count_aggr_df)
}
