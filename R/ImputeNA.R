
#' Title
#'
#' @param dt
#' @param val
#' @param na.imputation.method
#' @param hybrid.mar
#' @param hybrid.mnar
#' @param n_knn
#' @param tune.sigma
#' @param q
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' # import missForest, imputeLCMD
imputeNA <- function(dt,
                     val = "Abundance",
                     na.imputation.method = "none",
                     hybrid.mar = "KNN",
                     hybrid.mnar = "QRILC",
                     n_knn = 15, # knn imputation
                     tune.sigma = 0.5, # QRILC, MinProb imputation
                     q = 0.01, # MinDet, MinProb imputation
                     ...) { # ... for all parameters related to missForest

  if (na.imputation.method != "none") {

    # col.p <- intersect(c("SampleID", "Pool", "Condition", "Channel", "BioRep", "TechRep"), names(dt))
    col.p <- "Filename"
    row.p <- intersect(c("Protein", "Peptide", "Feature"), names(dt))
  mwk <- dt[, c(..row.p, ..col.p, ..val)]

  if ("Feature" %in% names(dt)) {
    res <- svDialogs::dlg_message(
      paste0(
        "Data is in Feature-level! The number of missing values can be very high! ",
        "By default the data will be rolled-up from Feature- to Peptide-level before ",
        "imputation. For having control on Feature-to-Peptide step, please check out the ",
        "vignette. \n\nDo you want to force imputation to be done on Feature level data? "
      ),
      "yesno"
    )$res
  } else {
    res <- "yes"
  }


  # long to wide/matrix
  if (res == "yes") {
    mwk <- MSqb2:::.l2w(
      dl = mwk, col.p = col.p, row.p = row.p, val = val, sp = "&.&.&",
      asMat = TRUE, Mat.rownm = row.p
    )
  } else {
    # !!!! TO DO: feature to peptide with median!!!! FARHAD: 17.07.2020
    # This must be changed using PSM2PPT function
    mwk <- MSqb2:::.l2w(
      dl = mwk, col.p = col.p, row.p = row.p, val = val, sp = "&.&.&",
      fun.agg = median, asMat = TRUE, Mat.rownm = row.p
    )
  }



  ## method QRILC
  if (toupper(na.imputation.method) == "QRILC") {
    imp.mat <- imputeLCMD::impute.QRILC(mwk, tune.sigma = tune.sigma)
    imp.mat <- imp.mat[[1]]
  }


  ## method KNN
  if (toupper(na.imputation.method) == "KNN") {
    imp.mat <- imputeLCMD::impute.wrapper.KNN(mwk, K = n_knn)
  }


  ## method MLE
  if (toupper(na.imputation.method) == "MLE") {
    imp.mat <- imputeLCMD::impute.wrapper.MLE(mwk)
  }


  ## method MinProb
  if (toupper(na.imputation.method) == "MINPROB") {
    imp.mat <- imputeLCMD::impute.MinProb(mwk, q = q, tune.sigma = tune.sigma)
  }


  ## method MinDet
  if (toupper(na.imputation.method) == "MINDET") {
    imp.mat <- imputeLCMD::impute.MinDet(mwk, q = q)
  }


  ## method zero
  if (toupper(na.imputation.method) == "ZERO") {
    imp.mat <- imputeLCMD::impute.ZERO(mwk)
  }


  ## method HYBRID
  if (toupper(na.imputation.method) == "HYBRID") {
    m.s <- model.Selector(mwk)
    imp.mat <- imputeLCMD::impute.MAR.MNAR(mwk, m.s,
      method.MAR = hybrid.mar,
      method.MNAR = hybrid.mnar
    )
  }


  ## method missForest
  if (toupper(na.imputation.method) %in% c("MISSFOREST", "RANDOMFOREST", "RANDOM FOREST", "RF")) {
    imp.mat <- missForest::missForest(mwk, verbose = TRUE, ...)
  }


  # wide to long
  imp.mat <- MSqb2:::.w2l(dw = imp.mat, col.p = col.p, row.p = row.p, val = val, sp = "&.&.&")
  return(imp.mat)
  } else return(dt)
}
