#' Median Normalization of Data
#'
#' This function performs median normalization on a dataset. Optionally, it can log2 transform the data before normalization
#' and perform normalization based on a subset of the data specified by the user.
#'
#' @param dat A `data.table` containing the data to be normalized. The table should include columns for the features, samples, and values to be normalized.
#' @param col.p A character string specifying the column name that contains the sample identifiers (e.g., `"SampleID"`).
#' @param row.p A character string specifying the column name that contains the feature identifiers (e.g., `"Protein"`).
#' @param val A character string specifying the column name that contains the values to be normalized (e.g., `"Abundance"`).
#' @param log2transform Logical. If `TRUE`, log2 transformation is applied to the data before normalization. Default is `TRUE`.
#' @param normalize.bySubset A character string specifying the column name to be used for subsetting the data before normalization. Default is `NULL`.
#' @param normalize.bySubset.matchTerm A character string specifying the term to match within the subset column for normalization. Default is `NULL`.
#' @param normalize.bySubset.matchTerm.exactMatch Logical. If `TRUE`, requires an exact match for the `normalize.bySubset.matchTerm`. Default is `FALSE`.
#'
#' @return A `data.table` containing the normalized data, with the same structure as the input `dat`.
#'
#' @details The `MEDnorm` function applies median normalization to the data, which adjusts the values in each sample such that
#' the median value is consistent across samples. This is particularly useful in MS or other omics data where
#' global shifts in signal intensities between samples need to be corrected.
#'
#' The function allows for:
#' \itemize{
#'   \item **Log2 Transformation**: Optional log2 transformation of the data prior to normalization.
#'   \item **Subset Normalization**: Normalization based on a subset of the data, defined by the `normalize.bySubset` parameter.
#' }
#' If `normalize.bySubset` is provided, the function calculates the median based on the specified subset, otherwise, it normalizes
#' the data using the entire dataset. The function also supports partial or exact matching for subsetting.
#'
#'
#' @examples
#' \dontrun{
#'   data <- data.table(
#'     SampleID = c("S1", "S2", "S3", "S4"),
#'     Protein = c("P1", "P2", "P1", "P2"),
#'     Abundance = c(1.2, 2.3, 3.1, 4.2)
#'   )
#'
#'   normalized_data <- MEDnorm(data, col.p = "SampleID", row.p = "Protein", val = "Abundance")
#' }
#'
#' @export
MEDnorm <- function(dat,
                    col.p,
                    row.p,
                    val,
                    log2transform = TRUE,
                    normalize.bySubset = NULL,
                    normalize.bySubset.matchTerm = NULL,
                    normalize.bySubset.matchTerm.exactMatch = FALSE) {

  mat <- dat[, unique(c(row.p, col.p, val, normalize.bySubset)), with = FALSE]

  # log2 transformation
  if (log2transform) mat[, (val) := log2(get(val))]


  # if normalization is to be done on a subset
  if (!is.null(normalize.bySubset)) {
    if (is.character(unlist(dat[, ..normalize.bySubset])) | is.factor(unlist(dat[, ..normalize.bySubset]))) { # charachter or factor
      if (is.null(normalize.bySubset.matchTerm)) {
        warning("normalize.bySubset.matchTerm is not defined. Normalization will be done using whole the dataset.")
        normalize.bySubset <- NULL
      } else {
        if (normalize.bySubset.matchTerm.exactMatch & !is.null(normalize.bySubset)) {
          medsmp <- mat[get(normalize.bySubset) == normalize.bySubset.matchTerm,
            .(medsmp = median(get(val), na.rm = TRUE)),
            by = col.p
          ]
        } else if (!is.null(normalize.bySubset)) {
          medsmp <- mat[get(normalize.bySubset) %like% normalize.bySubset.matchTerm,
            .(medsmp = median(get(val), na.rm = TRUE)),
            by = col.p
          ]
        }
      }
    } else if (!is.null(normalize.bySubset) && is.logical(unlist(mat[, ..normalize.bySubset]))) { # if normalize.bySubset is of class logical
      medsmp <- mat[which(get(normalize.bySubset)), .(medsmp = median(get(val), na.rm = TRUE)), by = col.p]
    } else {
      warning(paste(
        "normalize.bySubset column must either be of class character or logical. ",
        "Normalization will be done on the whole dataset."
      ))
    }
  }


  # if normalization is to be done on the whole dataset
  if (is.null(normalize.bySubset)) medsmp <- mat[, .(medsmp = median(get(val), na.rm = TRUE)), by = col.p]


  mat <- merge(mat, medsmp, by = col.p, all.x = TRUE)
  mat[, (val) := get(val) - medsmp + median(medsmp, na.rm = TRUE)]
  mat[, medsmp := NULL]

  if (ncol(dat) != ncol(mat)) {
    mat <- merge(dat[, -..val], mat, by = intersect(names(dat), c(row.p, col.p)))
  }

  MSqb2::char2fact(mat)
  return(mat)
}
