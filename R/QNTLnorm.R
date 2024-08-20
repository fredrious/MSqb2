
QNTLnorm <- function(dt,
                     col.p,
                     row.p,
                     val,
                     log2transform = TRUE) {
  mwk <- unique(dt[, c(..row.p, ..col.p, ..val)])
  # log2 transformation
  if (log2transform) mwk[, (val) := log2(get(val))]
  # long to wide/matrix
  mwk <- MSqb2:::.l2w(dl = mwk, col.p = col.p, row.p = row.p, val = val, asMat = TRUE, Mat.rownm = row.p)

  dimnm <- dimnames(mwk)
  mwk <- preprocessCore::normalize.quantiles(mwk)
  dimnames(mwk) <- dimnm

  # wide to long and merge
  mwk <- MSqb2:::.w2l(dw = mwk, col.p = col.p, row.p = row.p, val = val)
  if (ncol(dt) != ncol(mwk)) {
    mwk <- merge(dt[, -..val], mwk, by = intersect(names(dt), c(row.p, col.p)))
  }

  MSqb2::char2fact(mwk)
  return(mwk)
}
