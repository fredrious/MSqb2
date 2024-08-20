

#' Title
#'
#' @param dw
#' @param var
#' @param groupby
#' @param techrep
#' @param medianpolish.byPool
#' @param block
#' @param summarization.method
#' @param medianpolish.method
#' @param rm.total.missingness
#' @param Tables.path
#'
#' @return
#' @export
#'
#' @examples
Pep2Prot <- function(dw,
                     var = "Abundance",
                     groupby = "Protein",
                     techrep = techrep,
                     medianpolish.byPool = FALSE,
                     block = "Pool",
                     summarization.method = "tmp",
                     medianpolish.method = "column.eff",
                     rm.total.missingness = TRUE,
                     Tables.path = Tables.path) {
  
  # col.p <- intersect(c("SampleID", "Pool", "Condition", "Channel", techrep), names(dw))
  col.p <- "Filename"
  row.p <- intersect(c("Protein", "Peptide", "Feature"), names(dw))
  form <- paste0(paste(row.p, collapse = "+"), "~", paste(col.p, collapse = "+"))

  if (medianpolish.byPool) {
    dw.mp <- NULL
    for (im in unique(dw[, get(block)])) {
      sub.dw <- dw[get(block) == im, ]
      sub.dw.w <- dcast.data.table(sub.dw, formula = form, value.var = eval(var), sep = "&.&.&") # cast input data
      cols.wid <- names(sub.dw.w[, !..col.p]) # added columns after casting

      if (toupper(summarization.method) == "TMP") {
        dw.mp <- rbind(dw.mp, sub.dw.w[, .p2pMedpol(.SD, method = medianpolish.method), by = .(get(groupby)), .SDcols = cols.wid])
      }
      if (toupper(summarization.method) == "SUM") {
        dw.mp <- rbind(dw.mp, sub.dw.w[, .p2pSum(.SD), by = .(get(groupby)), .SDcols = cols.wid])
      }
    }
  } else {
    dw.w <- dcast.data.table(dw, formula = form, value.var = eval(var), sep = "&.&.&") # cast input data -- wide format
    cols.wid <- names(dw.w[, !..row.p]) # added columns after casting


    if (toupper(summarization.method) == "TMP") {
      dw.w <- dw.w[, .p2pMedpol(.SD, method = medianpolish.method), by = .(get(groupby)), .SDcols = cols.wid]
    }

    if (toupper(summarization.method) == "SUM") {
      dw.w <- dw.w[, .p2pSum(.SD), by = .(get(groupby)), .SDcols = cols.wid]
    }
    dw.mp <- dw.w
  }

  if ("get" %in% names(dw.mp)) setnames(dw.mp, old = "get", new = eval(groupby))

  if (length(col.p) > 1) {
    dw.mp[, eval(col.p) := tstrsplit(Channel, "&.&.&", fixed = TRUE)]
    if (!"Channel" %in% col.p) dw.mp[, Channel := NULL]
  } else {
    setnames(dw.mp, old = "Channel", new = eval(col.p))
  }
  MSqb2::char2fact(dw.mp)
  int2fact(dw.mp)



  ## remove rows with complete missingness
  if (rm.total.missingness) {
    dtw <- dcast(dw, formula = form, value.var = var)
    matw <- dtw[, -..col.p] %>% as.matrix(.)

    i.fullNA <- which(apply(matw, 1, function(x) sum(!is.na(x))) == 0)
    if (length(i.fullNA) != 0) dtw <- dtw[-i.fullNA, ]
    message(paste0(length(i.fullNA), " proteins were removed due to complete missingness across samples!"))
  }

  ## add number of peptides per protein
  if (groupby %in% names(dw)) {
    dw[, Peptide.count := uniqueN(Peptide), by = groupby]
    dw[, PSM.count := uniqueN(Feature), by = groupby]
    dw.mp <- merge(dw.mp, unique(dw[, c(..groupby, "Peptide.count", "PSM.count")]), by = groupby, all.x = TRUE)
  }


  pwb <- createWorkbook(title = "Protein Abundance Data")
  addWorksheet(pwb, "ProteinAbundance")
  writeDataTable(wb = pwb, x = dw.mp, sheet = "ProteinAbundance", colNames = TRUE)

  saveWorkbook(pwb, file = file.path(Tables.path, "ProteinAbundanceData.xlsx"), overwrite = TRUE)
  message("\nAn excel file with protein abundance data was saved in the following path:")
  message(file.path(Tables.path, "ProteinAbundanceData.xlsx"))

  MSqb2::char2fact(dw.mp)
  return(dw.mp)
}






.p2pSum <- function(dws, method) {
  clnm <- colnames(dws)
  dws <- as.matrix(dws)
  tmp <- apply(dws, 2, sum) # Abundance = overall median + column effect

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
