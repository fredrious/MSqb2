

#' Title
#'
#' @param feat.dt
#' @param metadata
#' @param protein.column
#' @param peptide.seq.column
#' @param peptide.annot.seq.column
#' @param charge.column
#' @param FileID.column
#' @param ion.score.column
#' @param quaninfo.column
#' @param isolation.interference.column
#' @param filename.column
#' @param abundance.columns
#' @param subsymDT
#'
#' @return
#' @export
#'
#' @examples
reformatRAW <- function(feat.dt = feat.dt,
                        metadata = metadata,
                        protein.column = "Master Protein Accessions",
                        peptide.seq.column = "Sequence",
                        peptide.annot.seq.column = "Annotated Sequence",
                        charge.column = "Charge",
                        FileID.column = "File ID",
                        ion.score.column = "Score",
                        quaninfo.column = "Quan Info",
                        isolation.interference.column = "Isolation Interference [%]",
                        desc.column = "Master Protein Descriptions",
                        contaminant.column = "Contaminant",
                        marked.as.column = "Marked as",
                        abundance.columns = "Abundance",
                        reference.channel = "reference",
                        ms.software = "PD",
                        subsymDT = data.table(
                          symout = c("\\~", "\\+", "\\:", "\\|", "\\-", "\\(", "\\)", "\\/", " "),
                          symin = c("_", ".p", ".", ".", "_", ".", ".", ".", ".")
                        )) {


  ## sbsetting the feat.dt to rquired columns
  parList <- c(
    Protein = protein.column,
    Peptide = peptide.seq.column,
    Peptide.annot = peptide.annot.seq.column,
    Charge = charge.column,
    FileID = FileID.column,
    Score = ion.score.column,
    Quaninfo = quaninfo.column,
    Interference = isolation.interference.column,
    Description = desc.column,
    MarkedAs = marked.as.column,
    Contaminant = contaminant.column
  )
  
  parList <- parList[parList %in% names(feat.dt)]
  if (!"Peptide" %in% names(parList) & "Peptide.annot" %in% names(parList)) {
    names(parList)[names(parList) == "Peptide.annot"] <- "Peptide"
  }
  
  
  ## guess method: TMT, LFQ
  if ("Channel" %in% names(metadata)) { ## if Channel is in metadata -> TMT
    MSmethod <- "TMT"
    ## Which columns contain intensity values
    abundance.columnss <- names(feat.dt)[names(feat.dt) %like% abundance.columns]
  } else { # LFQ
    MSmethod <- "LFQ"
    ## Which columns contain intensity values
    abundance.columnss <- "Intensity"
  }

  ## subset columns 
  feat.dt <- c(parList, abundance.columnss) %>% intersect(., names(feat.dt)) %>% feat.dt[, ., with = FALSE]


  ## cleaning names and parameters --> cleanDataPara returns 2 variable in the
  ## cln.nm.dt (clean-named data) and ref.name.tbl
  metadata %<>% cleanDataPara(
    ddt = .,
    subsymDT = subsymDT,
    complete.check = TRUE
  ) %>% .[["cln.nm.dt"]]

  feat.dt %<>% cleanDataPara(
    ddt = .,
    subsymDT = subsymDT,
    complete.check = FALSE
  )
  ref.name.tbl <- feat.dt[["ref.name.tbl"]]
  feat.dt <- feat.dt[["cln.nm.dt"]]


  ## record old, new and nicknames
  ref.name.tbl %<>%
    .[!original.names %like% abundance.columns, ] %>%
    merge(., as.data.table(parList, keep.rownames = "nicknames"),
      by.x = "original.names", by.y = "parList"
    )

  cat(
    "\n\nSpecial charachters in column names of the raw data and experimental",
    "design table were substituted using 'subsymDT' table.\n"
  )
  print(data.table(
    "Special character" = gsub("\\\\", "", subsymDT$symout),
    "Replaced by" = subsymDT$symin
  ))



  ## replace names with nicknames
  ref.name.tbl %<>% .[new.names %in% names(feat.dt), ]
  setnames(feat.dt,
    old = unlist(ref.name.tbl[, new.names]),
    new = unlist(ref.name.tbl[, nicknames])
  )


  ## extract gene names from description column
  if ("Description" %in% names(feat.dt)) {
    annotdata <- feat.dt[, .(Protein, Description)] %>% unique
    annotdata[, Genes := lapply(as.character(Description), \(x) MSqb2:::.psmDescGenes(x) ) %>% unlist ]
    feat.dt[, Description := NULL]
  } else annotdata <- NULL




  ## report proteomics strategy
  cat("\n\n---- Mass Spec method: ", MSmethod, " ----\n\n")
  if (MSmethod == "TMT") {
    ## extract channel names
    # channels <- unlist(ref.name.tbl[original.names %like% abundance.columns, new.names])
    channels.cols <- feat.dt %>% names() %>% setdiff(., ref.name.tbl$nicknames)
    channels <- gsub(paste0(abundance.columns, "|[[:punct:]]"), "", channels.cols)
    setnames(feat.dt, channels.cols, channels)
    
    cat("TMT channels:\n")
    print(data.frame(channels.cols, channels))
  }



  ## Basic statistics "Before Filtering"
  if (MSmethod == "TMT" | MSmethod == "LFQ") {
    ## count number of unique PSMs and define feature, Pools and fractions
    feat.dt <- unique(MSqb2::char2fact(feat.dt))
    if (!"Peptide.annot" %in% names(feat.dt)) {
      feat.dt[, Feature := factor(paste(Protein, Peptide, Charge, sep = "_"))]
    } else {
      feat.dt[, Feature := factor(paste(Protein, Peptide.annot, Charge, sep = "_"))]
    }
    ## Basic statistics "Before Filtering"
    dtcnt <- data.table(
      "Proteins (before Filtering)" = uniqueN(feat.dt$Protein),
      "Peptides (before Filtering)" = uniqueN(feat.dt$Peptide),
      "Features (before Filtering)" = uniqueN(feat.dt$Feature)
    )
  }


  ## wide to long format (TMT)
  if (MSmethod == "TMT") {
    id.vars <- setdiff(names(feat.dt), channels)
    feat.dt %<>% melt(.,
      value.name = "Intensity",
      variable.name = "Channel",
      id.vars = id.vars) %>% 
    MSqb2::char2fact(.)
  }



  ##  merge long-formatted data with metadata
  if (ms.software == "PD") {

    if (MSmethod == "TMT") {
      id.cols <- c("FileID", "Channel")
    } else if (MSmethod == "LFQ") {
      id.cols <- c("FileID")
    }
    
    # extract pool and fraction/FileID info
    extMixFrac <- decompFileID(feat.dt, MSmethod, id.cols)
    if (length(setdiff(names(extMixFrac), names(feat.dt))) > 0) {
      feat.dt <- merge(feat.dt, extMixFrac, by = id.cols)
    }
    
    int.nm <- intersect(names(feat.dt), names(metadata))

    if (length(int.nm) == 0) stop("Experiment design and MS data can not be merged!")

    if (MSmethod == "TMT") {
      feat.dt <- mergeMetadata(feat.dt, metadata, extMixFrac, int.nm)
      metadata <- feat.dt$metadata
      feat.dt <- feat.dt$msDT
    }

    if (MSmethod == "LFQ") {
      feat.dt <- merge(feat.dt, metadata, by = int.nm, all.y = TRUE)
    }
  }

  ## remove reference channel (!!! for now. Needs revision!)
  # if (Removereference.channel & MSmethod != "LFQ") {
  if (MSmethod != "LFQ" & !is.null(reference.channel)) {
      feat.dt <- feat.dt[!toupper(Channel) %in% toupper(reference.channel), ]
  }
  MSqb2::char2fact(feat.dt)

  ## table of all categorical variables
  if (exists("extMixFrac")) {
    metacols <- c(names(metadata), names(extMixFrac[, -"FileID"])) %>% unique()
  }
  metadata <- unique(feat.dt[, metacols, with = FALSE])
  setorderv(metadata, names(metadata))
  MSqb2::char2fact(metadata)

  ## add ms.software and MSmethod to parent environement
  list2env(list("ms.software" = ms.software, "MSmethod" = MSmethod), envir = .GlobalEnv)

  list("metadata" = metadata, "feat.dt" = feat.dt, "annotdata" = annotdata) %>%
  return(.)
}





## add Pool and Fraction !!! -> at the moment only for PD data
## This only works if column 'FileID' has the following format:
## FX.XX : X is Pool id and XX is the Fraction id. Exmpl: F2.11 means Pool 2, Fraction 11
## look at Proteome Discoverer User Guide, Software Version 2.2, section: Assigning the Order of Fractions
decompFileID <- function(msDT, MSmethod, id.cols) {
  
  dummy <- unique(msDT[, id.cols, with = FALSE]) %>% 
    .[, Filename := do.call(paste, c(.SD, sep = "_")), .SDcols = id.cols] %>% 
    .[order(as.integer(as.factor(Filename)))]
  
  
  if (MSmethod == "TMT") {
    flID <- "Pool"
    extMixFrac <-
      dummy %>%
      .[, c(flID, "Fraction") := tstrsplit(gsub("[[:alpha:]]", "", FileID), "[[:punct:]]")]
  } else if (MSmethod == "LFQ") {
    flID <- "SampleID"
    if (any(dummy$FileID %like% "[[:punct:]]")) {
      extMixFrac <- dummy %>%
        .[, c(flID, "Fraction") := tstrsplit(FileID, "[[:punct:]]")]
    } else { # no fractionation
      extMixFrac <- dummy %>%
        .[, Fraction := 1] %>%
        .[, SampleID := FileID]
    }
  }
  
  
  
  if (sum(is.na(extMixFrac$Fraction)) > 0) {
    message(
      "\n!!! The following FileIDs have no info on fractionation:",
      as.charachter(extMixFrac[is.na(Fraction), FileID]),
      "The Fraction IDs will be defined by integer sequencing along ordered
      spectrum files per 'Pool'.\n"
    )
    
    extMixFrac[is.na(Fraction), Fraction := paste0("f", seq_along(Filename)), by = "FileID"]
  }
  
  
  isFrac <- extMixFrac[, .(frc = uniqueN(Filename)), by = flID]
  if (all(isFrac$frc == 1)) {
    message("\n\n---No fractionation!")
    extMixFrac[, Fraction := NULL]
  } else {
    if (uniqueN(extMixFrac$Fraction) < 10) {
      extMixFrac[, Fraction := paste0("f", Fraction)]
    } else {
      extMixFrac[, Fraction := paste0("f", sprintf("%02s", Fraction))]
    }
    
    cat("\nTable of fractions:\n")
    extMixFrac[, -"Filename"] %>% unique() %>% print() 
  }
  
  if (MSmethod == "TMT") {
    if (uniqueN(extMixFrac$Pool) == 1) {
      extMixFrac[, Pool := NULL]
    } else extMixFrac[, Pool := paste0("Pool", Pool)]
  }
  
  return(MSqb2::char2fact(extMixFrac))
}
