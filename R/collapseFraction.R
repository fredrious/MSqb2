#'  Fraction combination
#'
#' @param dt
#' @param group
#' @param method
#' @param path
#'
#' @return
#' @export
#'
#' @examples
collapseFractions <- function(dt = dt,
                              group = "Pool",
                              method = "MAXset",
                              path = Tables.path) {
  
  if (!group %in% names(dt)) group <- NULL
  
  if ("Fraction" %in% names(dt)) {
    count.pre <- dt[, c(
      .("Protein count" = uniqueN(Protein)),
      .("Peptide count" = uniqueN(Peptide)),
      .("when" = "before Comb.")
    ),
    by = c("Fraction", group)
    ]


    by.col <- unique(c("Feature", group))
    ## if TMT use channel, if LFQ use condition
    if (!"Channel" %in% names(dt)) ch.par <- "Condition" else ch.par <- "Channel"

    if (!is.null(dt$Abundance)) {
      setnames(dt, "Abundance", "Intensity")
      message("\n---Fractions are being combined after normalization.\n")
      when <- "after"
      renameAbund2Int <- TRUE
    } else {
      message("\n---Fractions are being combined before normalization.\n")
      when <- "before"
      renameAbund2Int <- FALSE
    }

    ## unique measurements (Feature only measured in 1 Fraction per Pool)
    sharedPSM.Mix <- dt[, c(.(cnt = uniqueN(Fraction))), by = by.col] # table(sharedPSM.Mix$cnt)
    wch0 <- dt[sharedPSM.Mix[cnt == 1, ], on = by.col][, -"cnt"]




    ## non-unique measurements (Feature measured in more than 1 Fraction per Pool)
    wch <- dt[sharedPSM.Mix[cnt > 1, ], on = by.col][, -"cnt"]

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
    dt <- unique(rbindlist(list(wch0, wch)))
    MSqb2::char2fact(dt)



    ## plot count before and after fraction combination
    if (toupper(method) != "SUM") {
      count.post <- dt[, c(
        .("Protein count" = uniqueN(Protein)),
        .("Peptide count" = uniqueN(Peptide)),
        .("when" = "after Comb.")
      ), by = c("Fraction", group)]
      count.dt <- rbindlist(list(count.pre, count.post), use.names = TRUE)
      MSqb2::char2fact(count.dt)
      int2fact(count.dt)

      count.dt <- melt(count.dt,
        id.vars = c("Fraction", group, "when"),
        variable.name = "what", value.name = "count"
      )

      count.dt[, count := lapply(count, as.numeric)]
      count.dt$count <- unlist(count.dt$count)


      if (file.exists(file.path(path, "DataSummary.xlsx"))) {
        swb <- loadWorkbook(file.path(path, "DataSummary.xlsx"))
        if ("PrtPepCount_fraction" %in% names(swb)) removeWorksheet(swb, "PrtPepCount_fraction")
        addWorksheet(swb, "PrtPepCount_fraction")
        writeDataTable(swb, "PrtPepCount_fraction", count.dt, colNames = TRUE)
        saveWorkbook(swb, file = file.path(path, "DataSummary.xlsx"), overwrite = TRUE)
      }
    }


    if (when == "after") {
      setnames(dt, "Intensity", "Abundance")
    }

    outCl <- intersect(names(dt), c("Run", "Fraction", "Filename"))
    return(list("dt.fcomb" = dt[, -..outCl], "PrtPepCount_fraction" = count.dt))
  } else {
    message("No data Fractionation!")
    return(list("dt.fcomb" = dt, "PrtPepCount_fraction" = NULL))
  }
}
