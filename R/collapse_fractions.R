#' Collapse Fractions in Mass Spectrometry Data
#'
#' This function combines fractions in mass spectrometry data using various methods to summarize the intensities
#' of features across fractions. The function supports multiple methods for summarizing the fractions, including
#' maximum intensity, median intensity, sum of intensities, and tagging fractions.
#'
#' @param dat A `data.table` containing the MS data. The table should include columns such as `Feature`,
#'   `Fraction`, `Condition` if label-free, `Channel` if TMT, and `Intensity`.
#' @param group A character string indicating the grouping variable for the analysis (e.g., "Pool"). Default is "Pool".
#' @param method A character string specifying the method to use for combining fractions. Options include:
#'   \itemize{
#'     \item `"MAXset"`: Select the fraction with the minimum missingness and maximum intensity across conditions/channels.
#'     \item `"MAX"`: Select the fraction with the maximum intensity across conditions/channels, regardless of missingness.
#'     \item `"MED"` or `"MEDIAN"`: Use the median intensity across fractions.
#'     \item `"SUM"`: Sum the intensities across fractions.
#'     \item `"SINGLE"`: Use only peptides with a single fraction and remove others.
#'     \item `"MULTI"`, `"SEPERATE"`, `"TAGFRAC"`, `"TAGFRC"`: Tag features by fractions and pass to summarization.
#'   }
#'   Default is `"MAXset"`.
#' @param path A character string specifying the path to save the output summary tables. Default is `Tables.path`.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item `"dat.fcomb"`: A `data.table` containing the combined fractions data.
#'     \item `"PrtPepCount_fraction"`: A `data.table` summarizing the protein and peptide counts before and after fraction combination.
#'   }
#'
#' @details The function first identifies unique and non-unique measurements (features measured in more than one fraction per group).
#'   It then combines the fractions based on the specified method. The function can also generate a summary of the protein and
#'   peptide counts before and after fraction combination.
#'
#'
#' @export
collapse_fractions <- function(dat = dat,
                              group = "Pool",
                              method = "MAXset",
                              path = Tables.path) {

  if (!group %in% names(dat)) group <- NULL

  if ("Fraction" %in% names(dat)) {
    count.pre <- dat[, c(
      .("Protein count" = uniqueN(Protein)),
      .("Peptide count" = uniqueN(Peptide)),
      .("when" = "before Comb.")
    ),
    by = c("Fraction", group)
    ]


    by.col <- unique(c("Feature", group))
    ## if TMT use channel, if LFQ use condition
    if (!"Channel" %in% names(dat)) ch.par <- "Condition" else ch.par <- "Channel"

    if (!is.null(dat$Abundance)) {
      setnames(dat, "Abundance", "Intensity")
      message("\n---Fractions are being combined after normalization.\n")
      when <- "after"
      renameAbund2Int <- TRUE
    } else {
      message("\n---Fractions are being combined before normalization.\n")
      when <- "before"
      renameAbund2Int <- FALSE
    }

    ## unique measurements (Feature only measured in 1 Fraction per Pool)
    sharedPSM.Mix <- dat[, c(.(cnt = uniqueN(Fraction))), by = by.col] # table(sharedPSM.Mix$cnt)
    wch0 <- dat[sharedPSM.Mix[cnt == 1, ], on = by.col][, -"cnt"]




    ## non-unique measurements (Feature measured in more than 1 Fraction per Pool)
    wch <- dat[sharedPSM.Mix[cnt > 1, ], on = by.col][, -"cnt"]

    if (toupper(method) %in% c("MULTI", "SEPERATE", "TAGFRAC", "TAGFRC")) {
      ## features will be tagged by fractions and passed to summaritzation step
      wch[, Feature := paste0(Feature, "_", Fraction)]
    }

    if (toupper(method) %in% c("MAX", "MAXSET")) { # get the maximum fraction
      summary.wch <- wch[, c(
        .(n.nonNA = sum(!is.na(Intensity))),
        .(mean.Channels = mean(Intensity, na.rm = TRUE))
      ),
      by = c(by.col, "Fraction")
      ]

      setorderv(summary.wch, c("Feature", group, "n.nonNA", "mean.Channels"))

      wch.frc <- summary.wch[summary.wch[, .I[which.max(n.nonNA)], by = by.col]$V1]
      wch <- wch[wch.frc[, -c("n.nonNA", "mean.Channels")], on = c(by.col, "Fraction")]
    }


    if (toupper(method) %in% c("MED", "MEDIAN")) {
      wch <- wch[, Intensity := lapply(.SD, function(x) median(x, na.rm = TRUE)),
        .SDcols = "Intensity", by = c(by.col, ch.par)
      ]
    }


    if (toupper(method) == "SUM") {
      wch <- wch[, Intensity := lapply(.SD, sum), .SDcols = "Intensity", by = c(by.col, ch.par)]
      wch <- unique(wch[, Fraction := NULL])
      wch0 <- unique(wch0[, Fraction := NULL])
    }


    if (toupper(method) == "SINGLE") { # continue with Peptides with only 1 Fractions
      wch <- NULL
    }



    ## bind unique subset (wch0) and wch
    dat <- unique(rbindlist(list(wch0, wch)))
    char2fact(dat)



    ## plot count before and after fraction combination
    if (toupper(method) != "SUM") {
      count.post <- dat[, c(
        .("Protein count" = uniqueN(Protein)),
        .("Peptide count" = uniqueN(Peptide)),
        .("when" = "after Comb.")
      ), by = c("Fraction", group)]
      count.dat <- rbindlist(list(count.pre, count.post), use.names = TRUE)
      char2fact(count.dat)
      int2fact(count.dat)

      count.dat <- melt(count.dat,
        id.vars = c("Fraction", group, "when"),
        variable.name = "what", value.name = "count"
      )

      count.dat[, count := lapply(count, as.numeric)]
      count.dat$count <- unlist(count.dat$count)


      if (file.exists(file.path(path, "DataSummary.xlsx"))) {
        swb <- loadWorkbook(file.path(path, "DataSummary.xlsx"))
        if ("PrtPepCount_fraction" %in% names(swb)) removeWorksheet(swb, "PrtPepCount_fraction")
        addWorksheet(swb, "PrtPepCount_fraction")
        writeDataTable(swb, "PrtPepCount_fraction", count.dat, colNames = TRUE)
        saveWorkbook(swb, file = file.path(path, "DataSummary.xlsx"), overwrite = TRUE)
      }
    }


    if (when == "after") {
      setnames(dat, "Intensity", "Abundance")
    }

    outCl <- intersect(names(dat), c("Run", "Fraction", "Filename"))
    return(list("dat.fcomb" = dat[, -..outCl], "PrtPepCount_fraction" = count.dat))
  } else {
    message("No data Fractionation!")
    return(list("dat.fcomb" = dat, "PrtPepCount_fraction" = NULL))
  }
}
