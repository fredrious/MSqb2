
#' Title
#'
#' @param msDT
#' @param metadata
#' @param ms.software
#' @param save2file
#' @param Tables.path
#'
#' @return
#' @export
#'
#' @examples
mkQCtables <- function(msDT = msDT,
                       metadata = metadata,
                       ms.software = "PD",
                       save2file = TRUE,
                       Tables.path = Tables.path) {
  if (ms.software == "PD" & MSmethod %in% c("TMT", "LFQ")) {

    ## Fractions overlap: Features appear in multiple Fractions (per Pool)
    if (uniqueN(msDT$Fraction) > 1) {
      PSMsperFraction <- .FracOverlap(msDT, MSmethod)
    } else {
      PSMsperFraction <- NULL
    }
  }


  ##  merge long-formatted data with metadata
  if (MSmethod %in% c("TMT", "LFQ") & !is.null(metadata)) {
    smry.col <- intersect(names(msDT), unique(c(names(metadata), "Fraction")) )
    MissingPerFracChnl.w <- msDT[, c(
      .(cnt.NA = sum(is.na(Intensity))),
      .(perc.NA = round(100 * sum(is.na(Intensity)) / .N, 1)),
      .(cnt.Peptide = uniqueN(Peptide)),
      .(cnt.Protein = uniqueN(Protein))
    ),
    by = smry.col
    ]

    MissingPerFracChnl <- .missFracChnl(MissingPerFracChnl.w, metadata, smry.col)
  } else {
    MissingPerFracChnl <- NULL
  }



  # ## subset columns to those needed!
  # if (ms.software == "MQ") setnames(msDT, "Run", "Filename")
  # output.var <- c(
  #   "Protein", "Peptide", "Feature", "Pool", "Fraction", "Channel",
  #   "Condition", "BioRep", "TechRep", "Intensity", "Filename"
  # ) %>%
  #   c(., names(metadata)) %>%
  #   unique(.)
  # msDT <- msDT[, output.var[output.var %in% names(msDT)], with = FALSE]



  ## save all summary data (for qc) in file and export to global env.
  swb <- createWorkbook(title = "Data Summary")
  gl.objs <- c("PSMsperFraction", "MissingPerFracChnl", "extMixFrac")
  gl.objs <- as.list(intersect(gl.objs, ls()))

  exp.env <- new.env()
  sapply(gl.objs, function(x) {
    if (!is.null(get(x))) {
      addWorksheet(swb, x)
      writeDataTable(swb, x, as.data.frame(get(x)), colNames = TRUE)
      lout <- list(x = get(x)) %>% setNames(., x)
      list2env(lout, envir = exp.env)
    }
  })

  if (save2file) {
    saveWorkbook(swb, file = file.path(Tables.path, "DataSummary.xlsx"), overwrite = TRUE)
    message("\n\nA summary of the data (basic statistics, missingnes, ...) was saved in the following path:")
    message(file.path(Tables.path, "DataSummary.xlsx"))
  }

  # exp.env$filteredDT <- msDT
  return(as.list(exp.env))
}



## missingness per fraction and channel
.missFracChnl <- function(MissingPerFracChnl.w, metadata, smry.col) {
  MissingPerFracChnl <-
    sapply(MissingPerFracChnl.w, function(x) !is.factor(x)) %>%
    names(MissingPerFracChnl.w)[.] %>%
    MissingPerFracChnl.w[, (.) := lapply(.SD, function(x) as.numeric(x)), .SDcols = .] %>%
    melt.data.table(., id.vars = smry.col, variable.factor = TRUE) %>%
    MSqb2::char2fact(.)

  # MissingPerFracChnl <- merge(MissingPerFracChnl,
  #   metadata,
  #   by = intersect(
  #     names(MissingPerFracChnl),
  #     names(metadata)
  #   )
  # )

  MissingPerFracChnl <-
    intersect(c("Pool", "Channel", "Fraction"), smry.col) %>%
    setorderv(MissingPerFracChnl, .)
  return(MissingPerFracChnl)
}



## check fractions overlap
.FracOverlap <- function(msDT, MSmethod) {
  # distribution of PSMs across fractions
  if (MSmethod %in% c("TMT", "LFQ")) {
    if (!"Pool" %in% names(msDT)) msDT$Pool <- factor("dummypool")
    multFrac <- msDT[, .(cnt = uniqueN(Fraction)), by = list(Feature, Pool)]
    PSMsperFraction <- multFrac[, .N, by = list(cnt, Pool)] # as.data.table(table(multFrac$cnt))
    PSMsperFraction[, perc := round(N / sum(N), 4) * 100, by = Pool]
    PSMsperFraction[, txt := paste0("#", N, "\n", perc, "%")]
    MSqb2::char2fact(msDT)
    int2fact(msDT)
    PSMsperFraction$N <- as.numeric(PSMsperFraction$N)
    PSMsperFraction <- PSMsperFraction[order(Pool, cnt)]
  }
  if (msDT$Pool[1] == "dummypool") {
    msDT$Pool <- NULL
    PSMsperFraction$Pool <- NULL
  }
  return(PSMsperFraction)
}
