#' Generate Quality Control Summary Tables
#'
#' This function creates quality control (QC) summary tables from MS data and metadata. It calculates basic statistics,
#' missingness per fraction and channel, and checks for feature overlap across fractions. The results can be saved to an Excel file
#' and returned as a list for further analysis.
#'
#' @param dat A `data.table` containing the mass spectrometry data. This table should include columns such as `Feature`, `Protein`,
#'   `Peptide`, `Pool`, `Channel`, `Fraction`, and `Intensity`.
#' @param metadata A `data.table` containing metadata associated with the mass spectrometry data. This table may include columns such as
#'   `Pool`, `Channel`, and `Fraction`.
#' @param ms.software A character string specifying the mass spectrometry software used, either of `"PD"` and `"MQ"`. Default is `"PD"` (Proteome Discoverer).
#' @param save Logical. If `TRUE`, the QC tables are saved to an Excel file. Default is `TRUE`.
#' @param Tables.path A character string specifying the directory where the output file will be saved. Default is the value of `Tables.path`.
#'
#' @return A list containing the QC summary tables generated during the analysis, which may include:
#' \itemize{
#'   \item `"PSMsperFraction"`: A table showing the overlap of features across fractions.
#'   \item `"MissingPerFracChnl"`: A table summarizing the missingness of data per fraction and channel.
#'   \item `"extMixFrac"`: External mix fractions data (if applicable).
#' }
#'
#' @details The `create_qc_tables` function generates several summary tables that are essential for quality control in mass spectrometry-based
#' proteomics studies:
#'
#' \itemize{
#'   \item **Fractions Overlap**: Evaluates how features are distributed across multiple fractions, particularly useful for TMT or LFQ methods.
#'   \item **Missingness Analysis**: Provides a summary of missing data points per fraction and channel.
#' }
#'
#' The function merges the MS data with the provided metadata, computes the QC statistics, and optionally saves the
#' results in an Excel file for easy review and further analysis.
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeDataTable saveWorkbook
#' @importFrom magrittr %>%
#' @export
create_qc_tables <- function(dat = dat,
                       metadata = metadata,
                       ms.software = "PD",
                       save = TRUE,
                       Tables.path = Tables.path) {
  if (ms.software == "PD" & MSmethod %in% c("TMT", "LFQ")) {

    ## Fractions overlap: Features appear in multiple Fractions (per Pool)
    if (uniqueN(dat$Fraction) > 1) {
      PSMsperFraction <- .FracOverlap(dat, MSmethod)
    } else {
      PSMsperFraction <- NULL
    }
  }


  ##  merge long-formatted data with metadata
  if (MSmethod %in% c("TMT", "LFQ") & !is.null(metadata)) {
    smry.col <- intersect(names(dat), unique(c(names(metadata), "Fraction")) )
    MissingPerFracChnl.w <- dat[, c(
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

  if (save) {
    saveWorkbook(swb, file = file.path(Tables.path, "DataSummary.xlsx"), overwrite = TRUE)
    message("\n\nA summary of the data (basic statistics, missingnes, ...) was saved in the following path:")
    message(file.path(Tables.path, "DataSummary.xlsx"))
  }

  # exp.env$filteredDT <- dat
  return(as.list(exp.env))
}



## missingness per fraction and channel
.missFracChnl <- function(MissingPerFracChnl.w, metadata, smry.col) {
  MissingPerFracChnl <-
    sapply(MissingPerFracChnl.w, function(x) !is.factor(x)) %>%
    names(MissingPerFracChnl.w)[.] %>%
    MissingPerFracChnl.w[, (.) := lapply(.SD, function(x) as.numeric(x)), .SDcols = .] %>%
    melt.data.table(., id.vars = smry.col, variable.factor = TRUE) %>%
    char2fact(.)


  MissingPerFracChnl <-
    intersect(c("Pool", "Channel", "Fraction"), smry.col) %>%
    setorderv(MissingPerFracChnl, .)
  return(MissingPerFracChnl)
}



## check fractions overlap
.FracOverlap <- function(dat, MSmethod) {
  # distribution of PSMs across fractions
  if (MSmethod %in% c("TMT", "LFQ")) {
    if (!"Pool" %in% names(dat)) dat$Pool <- factor("dummypool")
    multFrac <- dat[, .(cnt = uniqueN(Fraction)), by = list(Feature, Pool)]
    PSMsperFraction <- multFrac[, .N, by = list(cnt, Pool)] # as.data.table(table(multFrac$cnt))
    PSMsperFraction[, perc := round(N / sum(N), 4) * 100, by = Pool]
    PSMsperFraction[, txt := paste0("#", N, "\n", perc, "%")]
    char2fact(dat)
    int2fact(dat)
    PSMsperFraction$N <- as.numeric(PSMsperFraction$N)
    PSMsperFraction <- PSMsperFraction[order(Pool, cnt)]
  }
  if (dat$Pool[1] == "dummypool") {
    dat$Pool <- NULL
    PSMsperFraction$Pool <- NULL
  }
  return(PSMsperFraction)
}
