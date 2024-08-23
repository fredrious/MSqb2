#' Quantile Normalization of Data
#'
#' This function performs quantile normalization on mass spectrometry data, with an option to log2 transform the data before normalization.
#'
#' @param dat A `data.table` containing the data to be normalized.
#' @param col.p Character vector. Column names used for grouping the data (e.g., sample or condition).
#' @param row.p Character vector. Row identifiers, typically the feature or protein IDs.
#' @param val Character. The column name representing the intensity or abundance values to be normalized.
#' @param log2transform Logical. Whether to log2 transform the data before normalization. Default is `TRUE`.
#'
#' @return A `data.table` with quantile normalized intensity values.
#'
#' @details
#' The function uses the `normalize.quantiles` function from the `preprocessCore` package to perform quantile normalization. This method makes the distribution of intensity values the same across samples. Log2 transformation is applied before normalization if `log2transform` is set to `TRUE`.
#'
#' @examples
#' \dontrun{
#'   normalized_data <- quantile_normalisation(dat = my_data, col.p = "SampleID", row.p = "ProteinID", val = "Intensity")
#' }
#'
#' @import preprocessCore
#' @export
quantile_normalisation <- function(dat,
                     col.p,
                     row.p,
                     val,
                     log2transform = TRUE) {
  mwk <- unique(dat[, c(..row.p, ..col.p, ..val)])
  # log2 transformation
  if (log2transform) mwk[, (val) := log2(get(val))]
  # long to wide/matrix
  mwk <- .l2w(dl = mwk, col.p = col.p, row.p = row.p, val = val, asMat = TRUE, Mat.rownm = row.p)

  dimnm <- dimnames(mwk)
  mwk <- preprocessCore::normalize.quantiles(mwk)
  dimnames(mwk) <- dimnm

  # wide to long and merge
  mwk <- .w2l(dw = mwk, col.p = col.p, row.p = row.p, val = val)
  if (ncol(dat) != ncol(mwk)) {
    mwk <- merge(dat[, -..val], mwk, by = intersect(names(dat), c(row.p, col.p)))
  }

  char2fact(mwk)
  return(mwk)
}
