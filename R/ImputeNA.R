#' Impute Missing Values in Data
#'
#' This function imputes missing values in a data set using various imputation methods. It supports a range of imputation techniques
#' including K-nearest neighbors (KNN), maximum likelihood estimation (MLE), quantile regression imputation for left-censored data (QRILC),
#' and hybrid methods combining different strategies for missing at random (MAR) and missing not at random (MNAR) data.
#'
#' @param dat A `data.table` containing the data to be imputed. The table should include columns corresponding to features, samples,
#'   and values to be imputed.
#' @param val A character string specifying the column name that contains the values to be imputed. Default is `"Abundance"`.
#' @param col.p A character string specifying the column name(s) that either contains unique sample names (e.g., `"Filename"`), or a combination of which can be used to uniqule address each sample (e.g., `c("Pool", "Channel")` for MS method TMT). Default is `"Filename"`.
#' @param row.p A character string specifying the column name that contains the feature identifiers. Default is `"Feature"`.
#' @param na.imputation.method A character string specifying the method to use for imputation. Options include:
#'   \itemize{
#'     \item `"V2"`: method v2 of the msImpute package
#'     \item `"V2-MNAR"`: method v2-mnar of the msImpute method specifically for MNAR data.
#'     \item `"QRILC"`: Quantile Regression Imputation of Left-Censored data.
#'     \item `"KNN"`: K-nearest neighbors.
#'     \item `"MLE"`: Maximum Likelihood Estimation.
#'     \item `"MINPROB"`: Impute using the minimal probability.
#'     \item `"MINDET"`: Impute using the minimal detection limit.
#'     \item `"ZERO"`: Impute with zeros.
#'     \item `"HYBRID"`: Hybrid method combining MAR and MNAR strategies.
#'   }
#'   Default is `"v2-mnar"`.
#' @param hybrid.mar A character string specifying the method to use for MAR imputation when using the hybrid approach. Default is `"KNN"`.
#' @param hybrid.mnar A character string specifying the method to use for MNAR imputation when using the hybrid approach. Default is `"QRILC"`.
#' @param n_knn An integer specifying the number of nearest neighbors to use for KNN imputation. Default is `15`.
#' @param tune.sigma A numeric value specifying the tuning parameter for the sigma in QRILC and MinProb methods. Default is `0.5`.
#' @param q A numeric value specifying the quantile for MinDet and MinProb methods. Default is `0.01`.
#' @param ... Additional arguments passed to the imputation functions.
#'
#' @return A `data.table` containing the original data with missing values imputed.
#'
#' @details The `imputeNA` function offers several methods for imputing missing data, depending on the nature of the missingness and the experimental design. For in-depth details please refere to vignette of the packages msImpute and imputeLCMD.
#'
#' \itemize{
#'   \item `"v2"`: This method, provided by the `msImpute` package, is designed for missing at random (MAR) data. It is particularly suitable for DIA experiments where missing values are assumed to be random.
#'   \item `"v2-mnar"`: This method extends `"v2"` to handle missing not at random (MNAR) data, making it suitable for cases where missing values are more likely for low-abundance features. This is often the case in DDA (like TMT) experiments, where the absence of data might be systematically related to the intensity of the signal.
#' }
#'
#'
#' @import imputeLCMD
#' @import msImpute
#' @export
imputeNA <- function(dat,
                     val = "Abundance",
                     col.p = "Filename",
                     row.p = "Feature",
                     na.imputation.method = "v2-mnar",
                     hybrid.mar = "KNN",
                     hybrid.mnar = "QRILC",
                     n_knn = 15, # for method knn
                     tune.sigma = 0.5, # for methods QRILC and MinProb
                     q = 0.01, # for methods MinDet, MinProb
                     ...) {
  # function body...
}

imputeNA <- function(dat,
                     val = "Abundance",
                     col.p = "Filename",
                     row.p = "Feature",
                     na.imputation.method = "v2-mnar",
                     hybrid.mar = "KNN",
                     hybrid.mnar = "QRILC",
                     n_knn = 15, # for method knn
                     tune.sigma = 0.5, # for methods QRILC and MinProb
                     q = 0.01, # for methods MinDet, MinProb
                     ...) {


  mat <- dat[, c(..row.p, ..col.p, ..val)] %>%
    MSqb2:::.l2w(dl = mat, col.p = col.p, row.p = row.p,
                 val = val, asMat = TRUE, Mat.rownm = "Filename")


  ## method v2 from msImpute
  if (toupper(na.imputation.method) == "V2") {
    imp.mat <- msImpute::msImpute(mat, method = "v2")
  }


  ## method v2 from msImpute
  if (toupper(na.imputation.method) == "V2-MNAR") {
    meta <- dat[, .(unique(c(col.p, "Condition"))), with = FALSE] %>% unique() %>%
      data.frame(., row.names = "Filename")
    group <- meta[colnames(mat), "Condition"] %>% as.vector() %>% as.character()

    imp.mat <- msImpute::msImpute(mat, method = "v2-mnar", group = group)
  }


  ## method QRILC
  if (toupper(na.imputation.method) == "QRILC") {
    imp.mat <- imputeLCMD::impute.QRILC(mat, tune.sigma = tune.sigma)
    imp.mat <- imp.mat[[1]]
  }


  ## method KNN
  if (toupper(na.imputation.method) == "KNN") {
    imp.mat <- imputeLCMD::impute.wrapper.KNN(mat, K = n_knn)
  }


  ## method MLE
  if (toupper(na.imputation.method) == "MLE") {
    imp.mat <- imputeLCMD::impute.wrapper.MLE(mat)
  }


  ## method MinProb
  if (toupper(na.imputation.method) == "MINPROB") {
    imp.mat <- imputeLCMD::impute.MinProb(mat, q = q, tune.sigma = tune.sigma)
  }


  ## method MinDet
  if (toupper(na.imputation.method) == "MINDET") {
    imp.mat <- imputeLCMD::impute.MinDet(mat, q = q)
  }


  ## method zero
  if (toupper(na.imputation.method) == "ZERO") {
    imp.mat <- imputeLCMD::impute.ZERO(mat)
  }


  ## method HYBRID
  if (toupper(na.imputation.method) == "HYBRID") {
    m.s <- model.Selector(mat)
    imp.mat <- imputeLCMD::impute.MAR.MNAR(mat, m.s,
      method.MAR = hybrid.mar,
      method.MNAR = hybrid.mnar
    )
  }


  # wide to long
  imp.mat <- MSqb2:::.w2l(dw = imp.mat, col.p = col.p, row.p = row.p, val = val)
  return(imp.mat)
}
