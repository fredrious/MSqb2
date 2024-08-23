#' Remove Batch Effects from Data
#'
#' This function removes batch effects from a dataset using the methods implimented in the packages `limma`, `sva`, or others. It returns the cleaned data with batch effects removed.
#'
#' @param dat A data.table containing the data to be corrected for batch effects.
#' @param metadt A data.table containing metadata corresponding to the data.
#' @param val A character string specifying the column name containing the values to be corrected. Default is "Abundance".
#' @param batch.corr.method A character string specifying the batch correction method to be used. Options are "limma" or "sva". Default is "limma".
#' @param row.p A character vector specifying the row identifiers in the data. Default is "auto".
#' @param col.p A character vector specifying the column identifiers in the data. Default is "auto".
#' @param model.formula A formula specifying the model for the batch correction.
#' @param sva.n.sv An integer specifying the number of surrogate variables to estimate if using the `sva` method. Default is `NULL`.
#'
#' @return A data.table with batch effects removed. If using `sva`, it returns a list containing the cleaned data and the `sva` object.
#' @importFrom sva sva
#' @importFrom magrittr %>% %<>%
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage of remove_batch_effect function
#' cleaned_data <- remove_batch_effect(dat, metadt, batch.corr.method = "limma")
#' }
remove_batch_effect <- function(dat,
                                metadt,
                                val = "Abundance",
                                batch.corr.method = "limma",
                                row.p = "auto",
                                col.p = "auto",
                                model.formula,
                                sva.n.sv = NULL) {

  # ComBat or removebatcheffects (limma package)
  # https://www.biostars.org/p/266507/


  # https://support.bioconductor.org/p/60879/
  # When calling removeBatchEffect, you should use the same design that you used for limma, but with
  # with batch effect term removed from the design. Then you would pass the batch effect factor as the
  # batch argument instead. So, if the design matrix that you used for limma was constructed as:
  #
  #   model.matrix(~Condition + Batch),
  #
  # then for removeBatchEffect, you would use design=model.matrix(~Condition), and batch=Batch.
  # In other words, you take the batch effect out of your model design and pass it as the batch
  # argument instead.

  sp <- "&.&.&" # dummy string for tstrsplt
  if (row.p[1] == "auto") row.p <- intersect(c("Protein", "Peptide", "Feature"), names(dat))
  if (col.p[1] == "auto") col.p <- names(metadt[, -"Fraction"])

  reftb <- dat[, ..col.p] %>% setorderv(., col.p) %>% unique() %>%
    .[, rr := do.call(paste, c(.SD, sep = sp)), .SDcols = col.p] %>% .[order(rr)]


  ## long to wide
  mt <- .l2w(
    dl = dat, col.p = col.p, row.p = row.p, val = val, sp = sp,
    asMat = TRUE, Mat.rownm = row.p
  ) %>% # order according to reftb
    .[, reftb$rr]



  ## method limma
  if (toupper(batch.corr.method) == "LIMMA") {

    dsgn <- model.matrix(as.formula(model.formula), data = as.data.frame(reftb))

    batch.vars <- model.formula %>% as.formula() %>% all.vars() %>% .[-1]
    if (length(batch.vars) > 0 ) batch1 <- batch.vars[1] else batch1 <- NULL
    if (length(batch.vars) > 1 ) batch2 <- batch.vars[2] else batch2 <- NULL

    # extract batch1 directly from the column names of the matrix
    if (!is.null(batch1)) {
      batch1 <- reftb[, ..batch1] %>% unlist() %>% unname() %>% as.character()
      # remove batch coluns from design matrix
      dsgn %<>% .[, !colnames(dsgn) %in% paste0(batch.vars[1], unique(batch1))]
    }
    if (!is.null(batch2)) {
      batch2 <- reftb[, ..batch2] %>% unlist() %>% unname() %>% as.character()
      # remove batch coluns from design matrix
      dsgn %<>% .[, !colnames(dsgn) %in% paste0(batch.vars[2], unique(batch2))]
    }

    cln.mat <- limma::removeBatchEffect(
      x = mt,
      batch = batch1,
      batch2 = batch2,
      design = dsgn
    )
  }



  ## method limma
  if (toupper(batch.corr.method) == "SVA") {

    # model formula must have only 1 variable. check for this!!!
    mod = model.matrix(as.formula(model.formula), data = metadt)
    mod0 = model.matrix(~1, data = metadt)

    # calculating the number of surrogate variables
    if (is.null(sva.n.sv)) sva.n.sv <- num.sv(mt, mod, method = "be")
    svobj = sva::sva(mt, mod, mod0, n.sv = sva.n.sv)

    # from: https://www.biostars.org/p/269093/
    Y <- t(mt)
    W <- svobj$sv
    alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
    svobj$corrected <- t(Y - W %*% alpha)
    cln.mat <- svobj$corrected
  }


  ## method limma
  # if (toupper(batch.corr.method) == "COMBAT") {
  # wk[, nNA := Reduce("+", lapply(.SD, function(x) is.na(x))), .SDcols = names(wk[, -..para])]
  # wk <- wk[nNA == 0]
  # wk <- wk[,!"nNA"]
  # }


  ## wide to long
  cln.dat <- .w2l(dw = cln.mat, col.p = col.p, row.p = row.p, val = val)

  ## merge with original data
  cln.dat %<>% merge(dat[, !..val], ., by = c(row.p, col.p))


  if (toupper(batch.corr.method) == "SVA") {
    list(cln.dat = cln.dat, svobj = svobj) %>% return()
  } else {
    return(cln.dat)
  }
}
