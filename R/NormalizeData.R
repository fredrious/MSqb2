

#' Title
#'
#' @param method
#' @param path
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' # import preprocessCore
NormalizeData <- function(method = "VSN", path = Tables.path, ...) {

  # args.norm <- list(...)
  list2env(list(...), environment())

  # col.p <- intersect(c("SampleID", "Pool", "Condition", "Channel", "TechRep", "Fraction"), names(dt))
  col.p <- "Filename"
  row.p <- intersect(c("Protein", "Peptide", "Feature"), names(dt))


  if (!is.null(normalize.bySubset)) {
    normalize.bySubset.matchTerm <- normalize.bySubset[2]
    normalize.bySubset.matchTerm.exactMatch <- ifelse(normalize.bySubset[3] == "exact", TRUE, FALSE)
    normalize.bySubset %<>% .[1]
  }


  if (toupper(method) == "QUANTILE" & !is.null(normalize.bySubset)) {
    warning(paste(
      "Normalising based on a data subset is only valid for methods VSN and MEDIAN. ",
      "Argument ubset will be ignored! Normalization will be done using the whole dataset."
    ))
  }

  switch(toupper(method),
    VSN = {
      msDTn <- VSNnorm(
        dt = dt, col.p = col.p, row.p = row.p, val = val, calib = calib,
        normalize.bySubset = normalize.bySubset,
        normalize.bySubset.matchTerm = normalize.bySubset.matchTerm,
        normalize.bySubset.matchTerm.exactMatch = normalize.bySubset.matchTerm.exactMatch
      )
    },
    MEDIAN = {
      msDTn <- MEDnorm(
        dt = dt, col.p = col.p, row.p = row.p, val = val,
        log2transform = TRUE,
        normalize.bySubset = normalize.bySubset,
        normalize.bySubset.matchTerm = normalize.bySubset.matchTerm,
        normalize.bySubset.matchTerm.exactMatch = normalize.bySubset.matchTerm.exactMatch
      )
    },
    QUANTILE = {
      msDTn <- QNTLnorm(
        dt = dt, col.p = col.p, row.p = row.p, val = val,
        log2transform = TRUE
      )
    }
  )


  # change intensity to abundance
  setnames(msDTn, old = val, new = "Abundance")


  # convert to wide to write in excel file
  dcast.form <- paste0(
    paste(row.p, collapse = "+"),
    "~",
    paste(names(msDTn[, -c(..row.p, "Abundance")]), collapse = "+")
  )

  nwb <- createWorkbook(title = "Normalised Data")

  addWorksheet(nwb, "Normalised_PeptideLevel")
  writeDataTable(
    wb = nwb, sheet = "Normalised_PeptideLevel", colNames = TRUE,
    x = dcast(msDTn, dcast.form, value.var = "Abundance")
  )

  addWorksheet(nwb, "NonNormalised_PeptideLevel")
  writeDataTable(
    wb = nwb, sheet = "NonNormalised_PeptideLevel", colNames = TRUE,
    x = dcast(dt, dcast.form, value.var = val)
  )

  saveWorkbook(nwb, file = file.path(path, "NormalisedData_PeptideLevel.xlsx"), overwrite = TRUE)

  message("\nAn excel file with normalised and non-normalised peptide-level data was saved in the following path:")
  message(file.path(path, "NormalisedData_PeptideLevel.xlsx"))

  return(msDTn)
}
