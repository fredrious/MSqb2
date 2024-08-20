
#' Title
#'
#' @param msd
#' @param metadata
#' @param filter.by.quaninfo
#' @param isolation.interference.cutoff
#' @param multi.features.method
#' @param filter.singleshot.proteins
#' @param min.intensity
#'
#' @return
#' @export
#'
#' @examples
filterFeatLevelData <- function(msd,
                                metadata,
                                filter.by.quaninfo = "auto",
                                isolation.interference.cutoff = 75,
                                multi.features.method = "TagPSM",
                                filter.singleshot.proteins = "BYPEPTIDE",
                                min.intensity = 0.01) {



  ## remove rows where protein is NA, SP, empty ...
  ##!!! depricated 27.2.2032 -> use "marked as" instead -> remove contaminants
  # msd %<>% .rmBadProt(.)


  ## Basic statistics "Before Filtering"
  np0 <- uniqueN(msd$Protein)
  dtcnt <- c(
    "Proteins (before Filtering)" = uniqueN(msd$Protein),
    "Peptides (before Filtering)" = uniqueN(msd$Peptide)
  )


  
  ## remove rows in which no accession code is avalable (missing)
  msd <- msd[!is.na(Protein)]
  
  
  ## remove contaminants
  if ("Contaminant" %in% names(msd)) {
    msd %<>% .[!(Contaminant)] %>%
      .[, Contaminant := NULL] %>% 
      .[!tolower(Protein) %like% "cont"]
  }
  if ("MarkedAs" %in% names(msd)) {
    msd %<>% .[!tolower(MarkedAs) %like% "contaminant"] %>%
      .[, MarkedAs := NULL] %>% 
      .[!tolower(Protein) %like% "cont"]
  }
  dtcnt["Proteins marked as contaminant"] <- np0 - uniqueN(msd$Protein)
  
  
  ## check single-shot proteins (before Peptide Filtering)
  msd %<>% .rmSnglProt(., filter.singleshot.proteins, when = "before", dtcnt)
  dtcnt <- c(dtcnt, msd[[2]])
  msd <- msd[[1]]


  ## remove zero and negative intensities
  msd %<>% .[Intensity <= min.intensity, Intensity := 1] %>% unique(.)


  cat("*** Peptide Filtering ***")

  ## exclude peptides based on info provided by quan.info
  if (filter.by.quaninfo != "none") {
    if ("Quaninfo" %in% names(msd)) {
      if (filter.by.quaninfo == "auto") {
        filter.by.quaninfo <- unique(msd[!is.na(Quaninfo), Quaninfo]) %>%
          toupper(.) %>%
          setdiff(., "UNIQUE")
      }
      badFeature <- unique(msd[toupper(Quaninfo) %in% filter.by.quaninfo, Feature])
      dtcnt["Miss.Quant/non-unique/No-Label Features"] = length(badFeature)

      cat("\n--- ", length(badFeature), " Features (~", round(length(badFeature) / uniqueN(msd$Feature), 2) * 100,
        "% of all Features) with the following labels will be removed.\n",
        paste(unique(filter.by.quaninfo), collapse = ", "),
        sep = ""
      )

      msd <- msd[!Feature %in% badFeature, ]
      msd$Quaninfo <- NULL
    }
  } else {
    dtcnt["Miss.Quant/non-unique/No-Label Features"] = "NA"
    badFeature <- c()
  }




  ## remove peptides by Isolation interference
  if (!is.null(isolation.interference.cutoff)) {
    if ("Interference" %in% names(msd)) {
      Isol.Interference <- list(
        "cutoff" = isolation.interference.cutoff,
        "Interference" = msd$Interference
      )
      n1 <- uniqueN(msd$Feature)

      cat("\n--- ", n1, " Features with Isolation Interference greater than ",
        isolation.interference.cutoff, "% will be removed!\n",
        sep = ""
      )


      msd <- msd[Interference < isolation.interference.cutoff, ]
      n2 <- uniqueN(msd$Feature)
      
      dtcnt[paste0(
        "Removed Feature due to high isolation interference (>",
        isolation.interference.cutoff, "%)"
      )] = n1 - n2
      
      msd$Interference <- NULL
    } else {
      message("!!! Unable to find Isolation Interference column in the data. Filter can not be applied!")
    }
  } else {
    dtcnt[paste0(
      "Removed redundant Feature due to high isolation interference (>",
      isolation.interference.cutoff, "%)"
    )] = "NA"
  }



  ## Dealing with multiple PSM measurements per feature and run.
  # msd %<>% MultiPSM(., method = multi.features.method, by.col = c(row.p, "Channel"))
  msd %<>% MultiPSM(.,
    method = multi.features.method, by.col =
      intersect(
        c("Protein", "Peptide", "Charge", "Filename", "Pool", "Channel"),
        names(msd)
      )
  )



  ## Filter peptides/features with complete missingness
  if (MSmethod == "TMT") {
    All_NA_Feat <- msd[, .(sm = sum(is.na(Intensity)), cn = .N), by = Feature] %>%
      .[sm == cn, Feature]
    msd %<>% .[!Feature %in% All_NA_Feat]

    cat("--- ", uniqueN(All_NA_Feat), " Features with total missingness will be removed!\n", sep = "")
  }



  ## Removing shared peptides between proteins (each peptide should belong only to 1 protein group)
    unqPpt <- msd[, .(n = uniqueN(Protein)), by = Peptide]

    if (nrow(unqPpt[n != 1]) != 0) {
      dtcnt["MissQuant/non-unique Peptides"] = nrow(unqPpt[n != 1]) + length(badFeature)
    }
    cat("--- ", nrow(unqPpt[n != 1]),
      " Peptides which are shared between multiple proteins (non-unique peptides) will be removed! ---",
      sep = ""
    )
    msd <- msd[Peptide %in% unqPpt[n == 1, Peptide], ]


  ## check single-shot proteins (after Peptide Filtering)
  msd %<>% .rmSnglProt(., filter.singleshot.proteins, when = "after", dtcnt)
  dtcnt <- c(dtcnt, msd[[2]])
  msd <- msd[[1]]


  
  ## continue with basic stats
  if (MSmethod %in% c("TMT", "LFQ")) {
    dtcnt["Features (after peptide Filtering)"] <- uniqueN(msd$Feature)
    # dtcnt["Removed proteins due to peptide Filtering"] <-
    #   dtcnt["Proteins (before Filtering)"] - 
    #   dtcnt["Single-shot Proteins (after peptide Filtering)"]
    dtcnt["Proteins (after Filtering)"] <- uniqueN(msd$Protein)
    BasicStats <- as.data.table(dtcnt, keep.rownames = TRUE)
    colnames(BasicStats) <- c("Parameter", "Count")
  }



  ## subset columns to those needed!
  if (ms.software == "MQ") setnames(msd, "Run", "Filename")
  output.var <- c(
    "Protein", "Peptide", "Feature", "Pool", "Fraction", "Channel",
    "Condition", "TechRep", "Intensity", "Filename", "Description"
  ) %>%
    c(., names(metadata)) %>%
    unique(.)
  msd <- msd[, output.var[output.var %in% names(msd)], with = FALSE]



  list(
    "msDT" = msd,
    "BasicStats" = BasicStats,
    "Isol.Interference" = Isol.Interference
  ) %>% return(.)
}




## find and remove SP Proteins!
## depricated 27.2.2032 -> use "marked as" instead -> remove contaminants
# .rmBadProt <- function(feat.dt) {
#   ## Removing SP proteins (rows with Master.Protein.Accessions == "SP")
#   feat.dt[, Protein := gsub(" ", "", Protein, fixed = TRUE)]
#   feat.dt[, Protein := toupper(Protein)]
#   # remove only SP protein groups (Master.Protein.Accessions)
#   # from rest of protein groups (Master.Protein.Accessions), remove SP term
#   feat.dt[, Protein := gsub("SP;", "", Protein, fixed = TRUE)]
#   feat.dt[, Protein := gsub(";SP", "", Protein, fixed = TRUE)]
#   feat.dt <- feat.dt[!Protein == "SP", ]
#   feat.dt <- feat.dt[!Protein == ";", ]
#   feat.dt[
#     startsWith(Protein, ";"),
#     Protein := lapply(Protein, function(x) substr(x, 2, nchar(x)))
#   ]
#   feat.dt[
#     endsWith(Protein, ";"),
#     Protein := lapply(Protein, function(x) substr(x, 1, nchar(x) - 1))
#   ]
#   return(feat.dt)
# }



## filter single shot proteins
.rmSnglProt <- function(dt, filter.singleshot.proteins, when = "before", dtcnt) {
  if (toupper(filter.singleshot.proteins) == "BYPEPTIDE") 
    bywhat <- ifelse("Peptide" %in% names(dt), "Peptide", "Peptide.annot")
  else 
    bywhat <- "Feature"

  unqPrt <- dt[, .(n = uniqueN(get(bywhat))), by = Protein]
  # dtcnt[, paste0("Single-shot Proteins (", when, " peptide Filtering)") := unqPrt[, sum(n == 1)]]
  info <- c(unqPrt[, sum(n == 1)])
  names(info) = paste0("Single-shot Proteins (", when, " peptide Filtering)")
  
  singleshot <- unqPrt[, 100 * round(sum(n == 1) / length(n), 4)]

  cat("\n\n*** single shot Proteins", toupper(when), "peptide Filtering ***")
  cat("\n--- ", unqPrt[, sum(n == 1)], " Proteins (",
    singleshot, "% of all proteins) have only 1 ", bywhat, " in the dataset ",
    sep = ""
  )
  cat(
    "\nfilter.singleshot.proteins = ", filter.singleshot.proteins,
    "\nsingle shot Proteins will be omitted from the dataset!\n\n"
  )
  dt[Protein %in% unqPrt[n != 1, Protein], ] %>%
    list(., info) %>%
    return(.)
}
