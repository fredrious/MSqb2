#' VSN Normalization
#'
#' This function applies Variance Stabilizing Normalization (VSN) to mass spectrometry data, with options to normalize based on subsets of the data or the entire dataset.
#'
#' @param dat A `data.table` containing the data to be normalized.
#' @param col.p Character vector. Column names used for grouping the data (e.g., sample or condition).
#' @param row.p Character vector. Row identifiers, typically the feature or protein IDs.
#' @param val Character. The column name representing the intensity or abundance values to be normalized.
#' @param calib Character. Calibration method to use, either `"affine"` (default) or `"none"`. If `"none"`, it assumes the data is already normalized.
#' @param normalize.bySubset Character or NULL. Column name used to define a subset for normalization. If `NULL`, the entire dataset is normalized.
#' @param normalize.bySubset.matchTerm Character or NULL. The term used to match entries within the subset defined by `normalize.bySubset`. If `NULL`, the whole dataset is used.
#' @param normalize.bySubset.matchTerm.exactMatch Logical. If `TRUE`, only exact matches of `normalize.bySubset.matchTerm` are used for subsetting. Default is `FALSE`.
#' @param sp Character. Separator used for pasting together row identifiers when converting data between wide and long formats. Default is `"_.._"`.
#'
#' @return A `data.table` with normalized intensity values.
#'
#' @details
#' This function uses the `vsn2` function from the `vsn` package to perform normalization. The normalization can be applied to the entire dataset or to a specific subset, which can be useful if different subsets require separate normalization.
#'
#' - **Subset Normalization:** If `normalize.bySubset` is specified, only the subset matching `normalize.bySubset.matchTerm` is used for normalization.
#' - **Calibration:** The `calib` parameter allows for affine calibration or no calibration (`"none"`).
#'
#' @examples
#' \dontrun{
#'   normalized_data <- vsn_normalisation(dat = my_data, col.p = "SampleID", row.p = "ProteinID", val = "Intensity")
#' }
#'
#' @import vsn
#' @import data.table
#' @export

vsn_normalisation <- function(dat,
                    col.p,
                    row.p,
                    val,
                    calib = "affine",
                    normalize.bySubset = NULL,
                    normalize.bySubset.matchTerm = NULL,
                    normalize.bySubset.matchTerm.exactMatch = FALSE,
                    sp = "_.._") {

  ## apply VSN
  # 1. calib can either be "affine" (default) or "none".
  #    "none" is usually being used when data is already calibrated (normalised, centered, ...).
  # 2. log2scale (from vsn manual): If TRUE, this will perform a global affine transform on the
  # data to put them on a similarscale as the original non-transformed data. Many users prefer
  # this. Fold-change estimates are not affected by this transform. In some situations, however,
  # it may be helpful to turn this off, e.g., when comparing independently normalized subsets of
  # the data.

  # use "none" when data is already normalized (i.e. by median). In this case, log2scale is
  # better set to FALSE!
  if (calib == "none") log2scale <- FALSE else log2scale <- TRUE


  # long to wide/matrix
  mwk <- dat[, unique(c(row.p, col.p, val, normalize.bySubset)), with = FALSE]
  mwk <- .l2w(dl = mwk, col.p = col.p, row.p = row.p, val = val, sp = sp)


  # if normalization is to be done on a subset
  # https://bioconductor.org/packages/release/bioc/vignettes/vsn/inst/doc/A-vsn.html
  if (!is.null(normalize.bySubset)) {
    if (is.character(unlist(dat[, ..normalize.bySubset])) | is.factor(unlist(dat[, ..normalize.bySubset]))) { # charachter or factor
      if (is.null(normalize.bySubset.matchTerm)) {
        warning("normalize.bySubset.matchTerm is not defined. Normalization will be done on whole dataset.")
        normalize.bySubset <- NULL
      } else {
        if (normalize.bySubset.matchTerm.exactMatch & !is.null(normalize.bySubset)) {
          subset.i <- which(mwk[, get(normalize.bySubset)] == normalize.bySubset.matchTerm)
        } else if (!is.null(normalize.bySubset)) {
          subset.i <- which(mwk[, get(normalize.bySubset)] %like% normalize.bySubset.matchTerm)
        }
        rownm <- mwk[, do.call(paste, c(.SD, sep = sp)), .SDcols = row.p]
        mwk <- as.matrix(mwk[, !..row.p], rownames = rownm)
        vsnfit <- vsn2(mwk[subset.i, ], calib = calib, verbose = FALSE)
      }
    } else if (!is.null(normalize.bySubset) && is.logical(unlist(mwk[, ..normalize.bySubset]))) { # if normalize.bySubset is of class logical
      subset.i <- which(mwk[, get(normalize.bySubset)])
      rownm <- mwk[, do.call(paste, c(.SD, sep = sp)), .SDcols = row.p]
      mwk <- as.matrix(mwk[, !..row.p], rownames = rownm)
      vsnfit <- vsn2(mwk[subset.i, ], calib = calib, verbose = FALSE)
    } else {
      warning(paste(
        "normalize.bySubset column must either be of class character or logical. ",
        "Normalization will be done using the whole dataset."
      ))
      normalize.bySubset <- NULL
    }
  }

  # if normalization is to be done on complete dataset
  if (is.null(normalize.bySubset)) {
    rownm <- mwk[, do.call(paste, c(.SD, sep = sp)), .SDcols = row.p]
    mwk <- as.matrix(mwk[, !..row.p], rownames = rownm)
    vsnfit <- vsn2(mwk, calib = calib, verbose = FALSE)
  }

  # normalised intensities
  pred <- vsn::predict(vsnfit, newdata = mwk, log2scale = log2scale)
  # VSNtest <- vsn::meanSdPlot(vsnfit, plot = FALSE)


  # wide to long and merge
  mwk <- .w2l(dw = pred, col.p = col.p, row.p = row.p, val = val, sp = sp)
  if (ncol(dat) != ncol(mwk)) {
    mwk <- merge(dat[, -..val], mwk, by = intersect(names(dat[, -..val]), c(row.p, col.p)))
  }

  char2fact(mwk)
  return(mwk)
}
