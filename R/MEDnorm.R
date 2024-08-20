

#' Title
#'
#' @param dt
#' @param col.p
#' @param row.p
#' @param val
#' @param log2transform
#' @param normalize.bySubset
#' @param normalize.bySubset.matchTerm
#' @param normalize.bySubset.matchTerm.exactMatch
#'
#' @return
#' @export
#'
#' @examples
MEDnorm <- function(dt,
                    col.p,
                    row.p,
                    val,
                    log2transform = TRUE,
                    normalize.bySubset = NULL,
                    normalize.bySubset.matchTerm = NULL,
                    normalize.bySubset.matchTerm.exactMatch = FALSE) {
  mwk <- dt[, unique(c(row.p, col.p, val, normalize.bySubset)), with = FALSE]
  # log2 transformation
  if (log2transform) mwk[, (val) := log2(get(val))]
  # long to wide/matrix
  # mwk <- MSqb2:::.l2w(dl = mwk, col.p = col.p, row.p = row.p, val = val)




  # if normalization is to be done on a subset
  if (!is.null(normalize.bySubset)) {
    if (is.character(unlist(dt[, ..normalize.bySubset])) | is.factor(unlist(dt[, ..normalize.bySubset]))) { # charachter or factor
      if (is.null(normalize.bySubset.matchTerm)) {
        warning("normalize.bySubset.matchTerm is not defined. Normalization will be done using whole the dataset.")
        normalize.bySubset <- NULL
      } else {
        if (normalize.bySubset.matchTerm.exactMatch & !is.null(normalize.bySubset)) {
          medsmp <- mwk[get(normalize.bySubset) == normalize.bySubset.matchTerm,
            .(medsmp = median(get(val), na.rm = TRUE)),
            by = col.p
          ]
        } else if (!is.null(normalize.bySubset)) {
          medsmp <- mwk[get(normalize.bySubset) %like% normalize.bySubset.matchTerm,
            .(medsmp = median(get(val), na.rm = TRUE)),
            by = col.p
          ]
        }
      }
    } else if (!is.null(normalize.bySubset) && is.logical(unlist(mwk[, ..normalize.bySubset]))) { # if normalize.bySubset is of class logical
      medsmp <- mwk[which(get(normalize.bySubset)), .(medsmp = median(get(val), na.rm = TRUE)), by = col.p]
    } else {
      warning(paste(
        "normalize.bySubset column must either be of class character or logical. ",
        "Normalization will be done on the whole dataset."
      ))
    }
  }


  # if normalization is to be done on the whole dataset
  if (is.null(normalize.bySubset)) medsmp <- mwk[, .(medsmp = median(get(val), na.rm = TRUE)), by = col.p]


  mwk <- merge(mwk, medsmp, by = col.p, all.x = TRUE)
  mwk[, (val) := get(val) - medsmp + median(medsmp, na.rm = TRUE)]
  mwk[, medsmp := NULL]

  if (ncol(dt) != ncol(mwk)) {
    mwk <- merge(dt[, -..val], mwk, by = intersect(names(dt), c(row.p, col.p)))
  }

  MSqb2::char2fact(mwk)
  return(mwk)
}
################################
