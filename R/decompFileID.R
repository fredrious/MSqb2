decompFileID2 <- function(dat, MSmethod) {
  
  # %%%%%%%%%%%%%%%%%%%%%%%%
  ## add Pool and Fraction !!! -> at the moment only for PD data
  ## This only works if column 'Run' has the following format:
  ## FX.XX : X is Pool id and XX is the Fraction id. Exmpl: F2.11 means Pool 2, Fraction 11
  ## look at Proteome Discoverer User Guide, Software Version 2.2, section: Assigning the Order of Fractions
  
  dummy <- unique(dat[, .(`File ID`, `Spectrum File`)]) %>% 
    MSqb2::order.num(., "File ID")
  
  if (MSmethod == "TMT") {
    flID <- "Pool"
    extMixFrac <-
      dummy %>%
      .[, c(flID, "Fraction") := tstrsplit(gsub("[[:alpha:]]", "", `File ID`), "[[:punct:]]")]
  } else if (MSmethod == "LFQ") {
    if (any(dummy$`File ID` %like% "[[:punct:]]")) {
      flID <- "SampleID"
      extMixFrac <- dummy %>%
        .[, c(flID, "Fraction") := tstrsplit(`File ID`, "[[:punct:]]")]
    } else { # no fractionation
      flID <- "File ID"
      extMixFrac <- dummy
    }
  }
  
  
  if (sum(is.na(extMixFrac$Fraction)) > 0) {
    message(
      "\n!!! The following `File ID`s have no info on fractionation:",
      as.charachter(extMixFrac[is.na(Fraction), `File ID`]),
      "The Fraction IDs will be defined by integer sequencing along ordered
      spectrum files per 'Pool'.\n"
    )
    extMixFrac[is.na(Fraction), Fraction := paste0("f", seq_along(`Spectrum File`)), by = "`File ID`"]
  }
  
  
  isFrac <- extMixFrac[, .(nfraction = uniqueN(`Spectrum File`)), by = flID]
  if (all(isFrac$nfraction == 1)) {
    message("\n\n---No fractionation!")
    if ("Fraction" %in% names(extMixFrac)) extMixFrac[, Fraction := NULL]
  } else {
    if (uniqueN(extMixFrac$Fraction) < 10) {
      extMixFrac[, Fraction := paste0("f", Fraction)]
    } else {
      extMixFrac[, Fraction := paste0("f", sprintf("%02s", Fraction))]
    }
  }
  
  if (MSmethod == "TMT") {
    if (uniqueN(extMixFrac$Pool) == 1) {
      if ("Pool" %in% names(extMixFrac))  extMixFrac[, Pool := NULL]
    } else extMixFrac[, Pool := paste0("Pool", Pool)]
  }
  
  extMixFrac <- MSqb2::order.num(extMixFrac, "File ID") %>% 
    MSqb2::char2fact() %>% 
    return()
}
