#' Summarize Peptide- or Feature-Level Data to Protein-Level Data
#'
#' This function aggregates peptide-level or Feature-Level data to protein-level data using various summarization methods.
#'
#' @param dat A data.table containing the peptide-level data.
#' @param var A character string specifying the column name that holds the values to be summarized (e.g., "Abundance").
#' @param groupby A character string specifying the column name used to group the data by proteins (e.g., "Protein").
#' @param summarization.method A character string specifying the summarization method to use. Options are `"tmp"` (median polish), `"sum"`, `"median"` and `"mean"`.f
#' @param medianpolish.byPool A logical indicating whether to perform median polish by pool.
#' @param block A character string specifying the column name used to define blocks (e.g., "Pool").
#' @param medianpolish.method A character string specifying the method for median polish summarization. Options are "column.eff" (default) or "residual".
#' @param rm.total.missingness A logical indicating whether to remove proteins with complete missingness across samples. Default is TRUE.
#' @param Tables.path A character string specifying the directory where the output Excel file with protein abundance data will be saved.
#'
#' @return A data.table containing the summarized protein-level data with additional columns for peptide and PSM counts.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # Example usage of the protein_aggregation function
#' protein_data <- protein_aggregation(
#'   dat = peptide_data,
#'   var = "Abundance",
#'   groupby = "Protein",
#'   medianpolish.byPool = FALSE,
#'   block = "Pool",
#'   summarization.method = "tmp",
#'   medianpolish.method = "column.eff",
#'   rm.total.missingness = TRUE,
#'   Tables.path = "path/to/save/tables"
#' )
#'}
protein_aggregation <- function(dat,
                                var = "Abundance",
                                groupby = "Protein",
                                summarization.method = "tmp",
                                medianpolish.byPool = FALSE,
                                block = "Pool",
                                medianpolish.method = "column.eff",
                                rm.total.missingness = TRUE,
                                Tables.path = Tables.path) {

  col.p <- "Filename"
  row.p <- intersect(c("Protein", "Peptide", "Feature"), names(dat))
  form <- paste0(paste(row.p, collapse = "+"), "~", paste(col.p, collapse = "+"))

  if (medianpolish.byPool) {
    dat.mp <- NULL
    for (im in unique(dat[, get(block)])) {
      sub.dat <- dat[get(block) == im, ]
      sub.dat.w <- dcast.data.table(sub.dat, formula = form, value.var = eval(var), sep = "&.&.&") # cast input data
      cols.wid <- names(sub.dat.w[, !..col.p]) # added columns after casting

      if (toupper(summarization.method) == "TMP") {
        dat.mp <- rbind(dat.mp, sub.dat.w[, .p2pMedpol(.SD, method = medianpolish.method), by = .(get(groupby)), .SDcols = cols.wid])
      } else {
        dat.mp <- rbind(dat.mp, sub.dat.w[, .p2pAgg(.SD, summarization.method), by = .(get(groupby)), .SDcols = cols.wid])
      }
    }
  } else {
    dat.w <- dcast.data.table(dat, formula = form, value.var = eval(var), sep = "&.&.&") # cast input data -- wide format
    cols.wid <- names(dat.w[, !..row.p]) # added columns after casting


    if (toupper(summarization.method) == "TMP") {
      dat.w <- dat.w[, .p2pMedpol(.SD, method = medianpolish.method), by = .(get(groupby)), .SDcols = cols.wid]
    } else {
      dat.w <- dat.w[, .p2pAgg(.SD, summarization.method), by = .(get(groupby)), .SDcols = cols.wid]
    }
    dat.mp <- dat.w
  }

  if ("get" %in% names(dat.mp)) setnames(dat.mp, old = "get", new = eval(groupby))

  if (length(col.p) > 1) {
    dat.mp[, eval(col.p) := tstrsplit(Channel, "&.&.&", fixed = TRUE)]
    if (!"Channel" %in% col.p) dat.mp[, Channel := NULL]
  } else {
    setnames(dat.mp, old = "Channel", new = eval(col.p))
  }
  char2fact(dat.mp)
  int2fact(dat.mp)



  ## remove rows with complete missingness
  if (rm.total.missingness) {
    dtw <- dcast(dat, formula = form, value.var = var)
    matw <- dtw[, -..col.p] %>% as.matrix(.)

    i.fullNA <- which(apply(matw, 1, function(x) sum(!is.na(x))) == 0)
    if (length(i.fullNA) != 0) dtw <- dtw[-i.fullNA, ]
    message(paste0(length(i.fullNA), " proteins were removed due to complete missingness across samples!"))
  }

  ## add number of peptides per protein
  if (groupby %in% names(dat)) {
    dat[, Peptide.count := uniqueN(Peptide), by = groupby]
    dat[, PSM.count := uniqueN(Feature), by = groupby]
    dat.mp <- merge(dat.mp, unique(dat[, c(..groupby, "Peptide.count", "PSM.count")]), by = groupby, all.x = TRUE)
  }


  pwb <- createWorkbook(title = "Protein Abundance Data")
  addWorksheet(pwb, "ProteinAbundance")
  writeDataTable(wb = pwb, x = dat.mp, sheet = "ProteinAbundance", colNames = TRUE)

  saveWorkbook(pwb, file = file.path(Tables.path, "ProteinAbundanceData.xlsx"), overwrite = TRUE)
  message("\nAn excel file with protein abundance data was saved in the following path:")
  message(file.path(Tables.path, "ProteinAbundanceData.xlsx"))

  char2fact(dat.mp)
  return(dat.mp)
}




.p2pAgg <- function(dws, method) {
  clnm <- colnames(dws)
  dws <- as.matrix(dws)
  tmp <- apply(X = dws, MARGIN = 2, FUN = method, na.rm = TRUE) # Abundance = overall median + column effect

  result <- data.table(Channel = clnm, Abundance = tmp)
  return(result)
}


.p2pMedpol <- function(dws, method) {
  clnm <- colnames(dws)
  dws <- as.matrix(dws)
  mp <- try(suppressWarnings(
    medpolish(dws,
              na.rm = TRUE,
              trace.iter = FALSE,
              maxiter = 50
    )
  ),
  silent = TRUE
  )

  if (toupper(method) %in% c("COLUMN.EFFECT", "COLUMN EFFECT", "COLUMN.EFF", "COLUMN")) {
    tmp <- mp$overall + mp$col # Abundance = overall median + column effect
  }

  if (toupper(method) %in% c("RES", "RESIDUAL")) {
    tmp <- mp$overall + apply(dws - mp$overall, 2, function(x) median(x, na.rm = TRUE)) # Abundance = overall median + column effect
  }

  result <- data.table(Channel = clnm, Abundance = tmp)
  # res <- mp$residual
  return(result)
}
