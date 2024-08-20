## cmplx is a list of multiple data.tables, each for one manually defined contrast
# complexContrast <- list(
#   cmplx.cont1 = data.table( ## case of even weighting, no double-delta
#     left    = c("cond1", "cond2"),
#     right   = c("cond3", "cond4")
#   ),
#   cmplx.cont2 = data.table( ## case of un-even weighting on the left
#     left    = c("cond3", "cond7", "cond6"),
#     right   = c("cond5"),
#     leftFC  = c(0.75, 0.3, 0.2)
#   ),
#   cmplx.cont3 = data.table( ## case of un-even weighting on the right
#     left    = c("cond3", "cond6"),
#     right   = c("cond5", "cond7"),
#     rightFC  = c(0.6, 0.4)
#   ),
#   deltadelta1 = data.table( ## case of double-delta with even weighting
#     left    = c("cond7", "cond6"),
#     right   = c("cond5", "cond3"),
#     leftFC  = c(1, -1),
#     rightFC = c(1, -1)
#   ),
#   deltadelta2 = data.table( ## case of double-delta with un-even weighting
#     left    = c("cond7", "cond6"),
#     right   = c("cond5", "cond3"),
#     leftFC  = c(0.75, -0.25),
#     rightFC = c(0.5, -0.5)
#   )
# )
## dds is the DESseq output: dds <- DESeq(dds)

makeComplexContrasts <- function(cmplx, contrastDT,
                              useComplexNames = FALSE,
                              prefix = NULL,
                              suffix = NULL,
                              export.as = "data.frame",
                              save2file = FALSE,
                              path = getwd()) {

  # contrastDT <- data.table(resNames = gsub(columns, "", resultsNames(dds) ))
  contrastDT <- unique(MSqb2::char2fact(contrastDT))
  dsgn.para <- names(contrastDT)

  for (idx in seq(length(cmplx))) {
    cmplxx <- as.data.table(cmplx[[idx]])
    MSqb2::char2fact(cmplxx)

    print(names(cmplx)[idx])
    print(cmplxx)

    nullFC <- setdiff(c("leftFC", "rightFC"), names(cmplxx))

    if (length(nullFC) == 1 & "rightFC" %in% nullFC) {
      message("rightFC is missing. Even-weighting on the right side of contrast formula.")
      cmplxx[, rightFC := 1 / uniqueN(cmplxx$right)]
      if (prod(sign(cmplxx$leftFC)) < 0) {
        message("leftFC has an odd number of negative factors (delta comparison!).\nFor delta-delta comparison both leftFC and rightFC must be provided.")
        stop()
      }
    }

    if (length(nullFC) == 1 & "leftFC" %in% nullFC) {
      message("leftFC is missing. Even weighting on the left side of contrast formula.")
      cmplxx[, leftFC := 1 / uniqueN(cmplxx$left)]
      if (prod(sign(cmplxx$rightFC)) < 0) {
        message("rightFC has an odd number of negative factors (delta comparison!).\nFor delta-delta comparison both leftFC and rightFC must be provided.")
        stop()
      } else if (prod(sign(cmplxx$rightFC)) > 0 & cmplxx$rightFC[1] > 0) {
        message("Provided rightFC are all posetive. It is assumed that the format in mind was -(rightFC.1 + rightFC.2 + ...)")
        message("All rightFC coeffients will therefore be converted to negative! Please check out the final contrast matrix.")
      }
    }

    if (length(nullFC) == 0) {
      if (any(c(prod(sign(cmplxx$rightFC)), prod(sign(cmplxx$leftFC))) < 0)) {
        message("leftFC and rightFC are both provided with negative factors. A delta-delta comparison!")
      }
    }

    if (length(nullFC) == 2) {
      cmplxx[, leftFC := 1 / uniqueN(cmplxx$left)]
      cmplxx[, rightFC := 1 / uniqueN(cmplxx$right)]
    }



    cmplxx[, ll := paste(leftFC, left, sep = "*")]
    cmplxx[, rr := paste(rightFC, right, sep = "*")]

    ## if even weighting, reomve 1* from string
    if (abs(prod(sign(cmplxx$leftFC))) == 1) cmplxx$ll <- gsub("1\\*", "", cmplxx$ll)
    if (abs(prod(sign(cmplxx$rightFC))) == 1) cmplxx$rr <- gsub("1\\*", "", cmplxx$rr)

    for (ir in 2:nrow(cmplxx)) {
      if (substring(cmplxx$ll[ir], 1, 1) != "-" & uniqueN(cmplxx$ll) != 1) cmplxx$ll[ir] <- paste0("+", cmplxx$ll[ir])
      if (substring(cmplxx$rr[ir], 1, 1) != "-" & uniqueN(cmplxx$rr) != 1) cmplxx$rr[ir] <- paste0("+", cmplxx$rr[ir])
    }

    contName <- paste(
      paste(unique(cmplxx$ll), collapse = ""),
      paste(unique(cmplxx$rr), collapse = ""),
      sep = " vs. "
    )

    cmplxx[, c("rr", "ll") := c(NULL, NULL)]

    if (!useComplexNames) contName <- names(cmplx)[idx] ## whethwe or not use long names for complex contrast

    ## after creating the names, wherever the rightFC set automatically, the coefficients must be negtive in cont. formula
    cmplxx$rightFC <- -1 * cmplxx$rightFC

    ## If there are common variables at the same position on each side, the final Coef.Fac. should be added...
    ## --> A-B vs. A-C must give C - B
    comPara <- intersect(cmplxx$left, cmplxx$right)
    if (length(comPara) > 0) {
      icx <- which(cmplxx == comPara, arr.ind = TRUE)
      CF <- cmplxx[icx[1, 1], icx[1, 2] + 2, with = FALSE] + cmplxx[icx[2, 1], icx[2, 2] + 2, with = FALSE]
      cmplxx[icx[1, 1], icx[1, 2] + 2] <- cmplxx[icx[2, 1], icx[2, 2] + 2] <- CF
    }

    contrastDT[, (contName) := 0]
    contrastDT <- merge(contrastDT, unique(cmplxx[, -c("left", "leftFC")]), by.x = names(contrastDT)[1], by.y = "right", all.x = TRUE)
    contrastDT <- contrastDT[!is.na(rightFC), (contName) := rightFC][, -"rightFC"]

    contrastDT <- merge(contrastDT, unique(cmplxx[, -c("right", "rightFC")]), by.x = names(contrastDT)[1], by.y = "left", all.x = TRUE)
    contrastDT <- contrastDT[!is.na(leftFC), (contName) := leftFC][, -"leftFC"]
  }

  message(
    "\n\n!!!WARNING...\nDOUBLE-CHECK THE COMPLEX CONTRASTS BELLOW. ALSO SAVED IN:\n",
    file.path(path, paste0(prefix, "complex_contrasts", suffix, ".xlsx")), "\n\n"
  )
  print(contrastDT)

  ## write contrasts into file
  if (save2file) {
    write.xlsx(
      as.data.table(contrastDT, keep.rownames = "Groups"),
      file.path(path, paste0(prefix, "complex_contrasts", suffix, ".xlsx")),
      overwrite = TRUE
    )
  }


  if (export.as == "data.frame") {
    rwnm <- as.character(unlist(contrastDT[, 1]))
    contrastDT <- as.data.frame(contrastDT[, -1])
    rownames(contrastDT) <- rwnm
    return(contrastDT)
  }
  if (export.as == "data.table") {
    return(contrastDT)
  }
}
