#' Collapse Multiple PSM Measurements per feature
#'
#' This function handles multiple PSM (Peptide Spectrum Match) measurements for a given feature, channel, and run by collapsing them into a single measurement according to the specified method.
#'
#' @param dat A `data.table` containing the PSM data.
#' @param method Character. The method used to collapse multiple PSM measurements. Options include `"IONSSCORE"`, `"MEAN"`, `"MEDIAN"`, `"MAX"`, and `"ALL"`. Default is `"IONSSCORE"`.
#' @param by Character vector. Columns used for grouping when collapsing PSM measurements.
#'
#' @return A `data.table` with collapsed PSM measurements according to the specified method.
#'
#' @details The function offers several methods for collapsing multiple PSM measurements of a given feature:
#' \itemize{
#'   \item \strong{IONSSCORE}: Selects the PSM with the highest ion score (ms method PD).
#'   \item \strong{MEAN}: Computes the mean intensity of all PSMs.
#'   \item \strong{MEDIAN}: Computes the median intensity of all PSMs.
#'   \item \strong{MAX}: Selects the maximum intensity among the PSMs.
#'   \item \strong{ALL}: Keeps all PSM measurements, tagging them for feature to protein summarization.
#' }
#'
#' @examples
#' \dontrun{
#'   collapsed_data <- collapse_psm(dat = my_data, method = "MEAN", by.col = c("Feature", "Channel", "Run"))
#' }
#'
#' @export
collapse_psm <- function(dat, method = "IONSSCORE", by) {
  if (toupper(method) %in% c("IONSSCORE", "IONSCORE", "ION.SCORE", "IONS.SCORE", "SCORE")) {
    if ("Score" %in% names(dat)) {
      dat <- dat[dat[, .I[which.max(Score)], by = by.col]$V1]
      dat$Score <- NULL
    } else {
      cat(
        "!!! Unable to find Ion.Score column in the data. ",
        "Filterring based on ion scores can not be applied!",
        "Instead, the redundant PSM measurements per Feature will ",
        "be tagged and passed to summarization step."
      )
      method <- "TAGPSM"
    }
  }
  if ("Score" %in% names(dat)) dat$Score <- NULL # not needed anymore


  if (toupper(method) == "MEAN") {
    dat[, Intensity := mean(Intensity, na.rm = TRUE), by = by.col] # Mean
  }

  if (toupper(method) %in% c("MED", "MEDIAN")) {
    dat[, Intensity := median(Intensity, na.rm = TRUE), by = by.col] # Median
  }

  if (toupper(method) == "MAX") {
    dat[, Intensity := max(Intensity, na.rm = TRUE), by = by.col] # max
  }

  if (toupper(method) %in% c("ALL")) {
    ## all spectrum level measurements (PSM in terms of PD) will enter MedianPlish procedure as an individual feature
    dat[, nPSM := length(Intensity), by = by.col]
    dat[nPSM > 1, Feature.ex := seq_along(Intensity), by = by.col]
    dat[nPSM > 1, Feature := paste0(Feature, "_PSM", Feature.ex)]
    dat[, nPSM := NULL]
    dat[, Feature.ex := NULL]
  }

  return(unique(char2fact(dat)))
}
