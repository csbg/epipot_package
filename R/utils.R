#' Remove special characters from a string.
#'
#' @param s The string to be cleaned.
#'
#' @return The string with parentheses/commas removed and all other
#'   non-alphanumeric characters replaced by underscores
clean_string <- function(s) {
  s <- gsub("[\\(\\),]", "", s)
  gsub("[^[:alnum:]]", "_", s)
}


#' Scale the rows of a matrix.
#'
#' @param m A numeric matrix
#'
#' @return A matrix whose rows have been scaled.
#' @export
#'
#' @examples
#' scale_rows(diag(4))
scale_rows <- function(m) {
  t(scale(t(m)))
}


#' Convert an annotated GRanges object to a data frame.
#'
#' This functions derives from `annoGR2DF()` in the Repitools package.
#'
#' @param gr A Granges object
#'
#' @returns A data frame describing genomic ranges.
#' @export
granges_to_dataframe <- function(gr) {
  df <- as.data.frame(gr)
  colnames(df)[1] <- "chr"
  if ('*' %in% df$strand)
    df <- df[, -match("strand", colnames(df))]
  df
}
