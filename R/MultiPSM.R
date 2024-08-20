
################################
#### --------------- mult. measurements per PSM
## Dealing with multiple PSM measurements perfeature, channel and run.

MultiPSM <- function(dt, method = "none", by.col) {
  if (toupper(method) %in% c("IONSSCORE", "IONSCORE", "ION.SCORE", "IONS.SCORE", "SCORE")) {
    if ("Score" %in% names(dt)) {
      dt <- dt[dt[, .I[which.max(Score)], by = by.col]$V1]
      dt$Score <- NULL
    } else {
      cat(
        "!!! Unable to find Ion.Score column in the data. ",
        "Filterring based on ion scores can not be applied!",
        "Instead, the redundant PSM measurements per Feature will ",
        "be tagged and passed to summarization step."
      )
      method <- "TAGPSM"
    }
  }
  if ("Score" %in% names(dt)) dt$Score <- NULL # not needed anymore


  if (toupper(method) == "MEAN") {
    dt[, Intensity := mean(Intensity, na.rm = TRUE), by = by.col] # Mean
  }

  if (toupper(method) %in% c("MED", "MEDIAN")) {
    dt[, Intensity := median(Intensity, na.rm = TRUE), by = by.col] # Median
  }

  if (toupper(method) == "MAX") {
    dt[, Intensity := max(Intensity, na.rm = TRUE), by = by.col] # max
  }

  if (toupper(method) %in% c("ALL")) {
    ## all spectrum level measurements (PSM in terms of PD) will enter MedianPlish procedure as an individual feature
    dt[, nPSM := length(Intensity), by = by.col]
    dt[nPSM > 1, Feature.ex := seq_along(Intensity), by = by.col]
    dt[nPSM > 1, Feature := paste0(Feature, "_PSM", Feature.ex)]
    dt[, nPSM := NULL]
    dt[, Feature.ex := NULL]
  }


  return(unique(MSqb2::char2fact(dt)))
}
################################
