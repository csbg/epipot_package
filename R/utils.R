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
#' @return A matrix whose rows have been scaled
#' @export
#'
#' @examples
scale_rows <- function(m) {
  t(scale(t(m)))
}
