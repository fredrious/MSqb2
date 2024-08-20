
VSNnorm <- function(dt,
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
  mwk <- dt[, unique(c(row.p, col.p, val, normalize.bySubset)), with = FALSE]
  mwk <- MSqb2:::.l2w(dl = mwk, col.p = col.p, row.p = row.p, val = val, sp = sp)


  # if normalization is to be done on a subset
  # https://bioconductor.org/packages/release/bioc/vignettes/vsn/inst/doc/A-vsn.html
  if (!is.null(normalize.bySubset)) {
    if (is.character(unlist(dt[, ..normalize.bySubset])) | is.factor(unlist(dt[, ..normalize.bySubset]))) { # charachter or factor
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
  mwk <- MSqb2:::.w2l(dw = pred, col.p = col.p, row.p = row.p, val = val, sp = sp)
  if (ncol(dt) != ncol(mwk)) {
    mwk <- merge(dt[, -..val], mwk, by = intersect(names(dt[, -..val]), c(row.p, col.p)))
  }

  MSqb2::char2fact(mwk)
  return(mwk)
}
