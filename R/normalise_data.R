#' Normalize MS Data
#'
#' This function normalizes mass spectrometry data using one of three available methods: Variance Stabilizing Normalization (VSN), median normalization, or quantile normalization. The function can also handle normalization by specific data subsets.
#'
#' @param method Character. The normalization method to use. Options are `"vsn"`, `"median"`, or `"quantile"`. Default is `"VSN"`.
#' @param path Character. The path where the output Excel file with normalized and non-normalized data will be saved. Default is `Tables.path`.
#' @param ... Additional arguments passed to the normalization functions, including the column and row identifiers (`col.p`, `row.p`), the intensity column (`val`), and other specific parameters for each normalization method.
#'
#' @return A `data.table` containing the normalized data, with the intensity values replaced by the normalized abundance values.
#'
#' @details
#' - **vsn**: Variance Stabilizing Normalization is applied, with an option to normalize by a subset of the data.
#' - **median**: Median normalization can be applied, with log2 transformation by default.
#' - **quantile**: Quantile normalization is performed, transforming the intensity distributions of all samples to be identical.
#'
#' The function also generates an Excel file containing both the normalized and non-normalized peptide-level data.
#'
#' @examples
#' \dontrun{
#' normalized_data <- normalise_data(
#'   method = "VSN",
#'   dat = my_data,
#'   val = "Intensity",
#'   path = "output_directory"
#' )
#' }
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeDataTable saveWorkbook
#' @importFrom magrittr %<>%
#' @export

normalise_data <- function(method = "VSN", path = Tables.path, ...) {

  # args.norm <- list(...)
  list2env(list(...), environment())

  # col.p <- intersect(c("SampleID", "Pool", "Condition", "Channel", "TechRep", "Fraction"), names(dat))
  col.p <- "Filename"
  row.p <- intersect(c("Protein", "Peptide", "Feature"), names(dat))


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
           msDTn <- vsn_normalisation(
             dat = dat, col.p = col.p, row.p = row.p, val = val, calib = calib,
             normalize.bySubset = normalize.bySubset,
             normalize.bySubset.matchTerm = normalize.bySubset.matchTerm,
             normalize.bySubset.matchTerm.exactMatch = normalize.bySubset.matchTerm.exactMatch
           )
         },
         MEDIAN = {
           msDTn <- median_normalisation(
             dat = dat, col.p = col.p, row.p = row.p, val = val,
             log2transform = TRUE,
             normalize.bySubset = normalize.bySubset,
             normalize.bySubset.matchTerm = normalize.bySubset.matchTerm,
             normalize.bySubset.matchTerm.exactMatch = normalize.bySubset.matchTerm.exactMatch
           )
         },
         QUANTILE = {
           msDTn <- quantile_normalisation(
             dat = dat, col.p = col.p, row.p = row.p, val = val,
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
    x = dcast(dat, dcast.form, value.var = val)
  )

  saveWorkbook(nwb, file = file.path(path, "NormalisedData_PeptideLevel.xlsx"), overwrite = TRUE)

  message("\nAn excel file with normalised and non-normalised peptide-level data was saved in the following path:")
  message(file.path(path, "NormalisedData_PeptideLevel.xlsx"))

  return(msDTn)
}
