

#' Title
#'
#' @param dat
#' @param feature.annotation.source
#' @param uniprot.input
#' @param uniprot.output
#' @param uniprot.extra
#' @param uniprot.KEY
#' @param uniprot.species.name
#' @param uniprot.taxId
#' @param n.call
#' @param uniprot.file
#' @param annotation.file
#' @param PD.Desc
#' @param PD.Acc
#' @param shorten.gene.names
#' @param split.char
#' @param save.annotation
#' @param save.annotation.path
#' @param Data.path
#'
#' @return
#' @export
#'
#' @examples
#'
#' # import UniProt.ws, vsvDialogs
callAnnotations <- function(dat,
                           feature.annotation.source, # "uniprot.file" , "annotation.file"
                           # if source "uniprot.web":
                           uniprot.input = "Protein",
                           uniprot.output = "Genes",
                           uniprot.extra = NULL, # more parameters can be called from uniprot
                           uniprot.KEY = "UNIPROTKB",
                           uniprot.species.name = NULL,
                           uniprot.taxId = NULL,
                           n.call = 99,
                           uniprot.file = "AnnotationsTable.xlsx",
                           shorten.gene.names = TRUE, # shorten the gene and protein names?
                           split.char = ";",
                           save.annotation = TRUE,
                           save.annotation.path = Tables.path,
                           Data.path = Data.path) {
  

  if (file.exists(feature.annotation.source)) src <- "annotation.file" else src <- "nofile"

  MasterPrtAll <- unique(dat[, ..uniprot.input])
  prot.cnt <- as.data.table(table(MasterPrtAll[, lengths(strsplit(as.character(get(uniprot.input)), ";"))]))


  if (sum(as.integer(prot.cnt$V1)) > 1) { # if exists inputs with more than 1 accession code (usually seperated by ';' )
    MasterPrtAll[, cntProt := lengths(strsplit(as.character(get(uniprot.input)), split.char))]

    ## subdata with only 1 accession code per input
    protDT <- MasterPrtAll[cntProt == 1][, paste0(uniprot.input, "1") := get(uniprot.input)]

    ## subdata with np number of accession codes per input
    nAcc <- unique(MasterPrtAll$cntProt)
    for (np in nAcc[nAcc > 1]) {
      tmpdt <- MasterPrtAll[cntProt == np]
      tmpdt[, paste0("Protein", seq(np)) := tstrsplit(get(uniprot.input), split.char, fixed = TRUE)]
      protDT <- rbind(protDT, tmpdt, fill = TRUE)
      rm(tmpdt)
    }

    cols <- paste0(uniprot.input, seq(max(nAcc)))
    protDT. <- protDT[, ..cols]
    MSqb2::char2fact(protDT.)

    # calling gene annotations from uniprot API
    allACCs <- data.table(Protein = unique(unlist(protDT.)[!is.na(unlist(protDT.))]))
  } else {
    allACCs <- unique(MasterPrtAll[, ..uniprot.input])
  }


  if (src == "psm.file") {
    gnDT <- 
      dat[, Genes := lapply(Descriptions, .psmDescGenes) %>% unlist ]
    
        # desc.dt[, c(Protein, Description)] %>%
        #   .[, tmp := tstrsplit(Description, "GN=")[2]] %>%
        #   .[, Genes := tstrsplit(tmp, "PE=")[1]] %>%
        #   .[, Genes := gsub(" ", "", Genes)] %>%
        #   .[, c("tmp", "Description") := NULL]
  } else {
    
    if (feature.annotation.source == "uniprot.file") {
      message(
        "When 'feature.annotation.source' set to 'uniprot.file', the last already-saved 'AnnotationsTable.xlsx' ",
        "file from the last call to uniprot website will be used for gene annotation."
      )
      
      fl <- list.files(
        path = save.annotation.path,
        pattern = "AnnotationsTable.xlsx",
        include.dirs = TRUE, ignore.case = TRUE,
        full.names = FALSE, recursive = TRUE
      )
      (fl <- file.info(fl, extra_cols = FALSE))
      fl <- file.path(rownames(fl[fl$ctime == max(fl$ctime), ]))
      
      if (length(fl) != 0) {
        message("selected annotation file:\n", fl)
        gnDT <- as.data.table(read_excel(fl))
      } else {
        message("Local file for gene annotation not found. Gene names wil be called directly from uniprot!")
        feature.annotation.source == "uniprot.web"
      }
      ## when reading annotations from a local file, there is no need to save the annotation file again!
      save.annotation <- FALSE
    }
    
    
    if (feature.annotation.source == "uniprot.web") {
      if (is.null(uniprot.taxId) & is.null(uniprot.species.name)) {
        stop(
          "Either of uniprot taxonomy ID or full uniprot species name",
          "must be provided!"
        )
      } else {
        if (is.null(uniprot.taxId) & !is.null(uniprot.species.name)) {
          taxonID <- availableUniprotSpecies(pattern = uniprot.species.name, n = Inf)
          
          if (nrow(taxonID) > 1) {
            message("Multiple matches with selected Species name! Please enter taxonomy ID from the list below:")
            print(taxonID)
            uniprot.taxId <- unlist(dlg_form(list("Taxonomy ID" = "Example: 10090"), "Enter Taxonomy ID:")$res)
            # uniprot.taxId <- dlg_input("Multiple matches with selected Species name! Please enter taxonomy ID from the list:")$res
            if (!length(uniprot.taxId)) {
              message("No taxonomy ID entered! Gene annotation will be skipped!\n")
            } else if (!uniprot.taxId %in% taxonID$"taxon ID") {
              message("Enetered taxonomy ID does not match any taxonomy ID in the list! Please try again...")
              # uniprot.taxId <- dlg_input("Please enter taxonomy ID:")$res
              uniprot.taxId <- unlist(dlg_form(list("Taxonomy ID" = "Example: 10090"), "Enter Taxonomy ID:")$res)
            }
          } else {
            uniprot.taxId <- taxonID$"taxon ID"
          }
        }
        
        if (!is.null(uniprot.taxId) & is.null(uniprot.species.name)) {
          uniprot.species.name <- lookupUniprotSpeciesFromTaxId(uniprot.taxId)
        }
        
        message(paste("\n\n** uniProt Species:", uniprot.species.name, "**"))
        message(paste("** uniProt taxonomy ID:", uniprot.taxId, "** \n\n"))
      }
      
      message(paste0("calling annotation for ", uniqueN(allACCs), " accession codes ... "))
      
      uniprot.taxId <- as.numeric(uniprot.taxId)
      uniprot.AnnotDB <- UniProt.ws(uniprot.taxId)
      
      .calluniprot <- function(input, gnDT) {
        gnDT <- rbindlist(
          list(
            gnDT,
            UniProt.ws::select(uniprot.AnnotDB,
                               keys = input,
                               taxId = uniprot.taxId,
                               columns = toupper(c(uniprot.output, uniprot.extra)),
                               keytype = uniprot.KEY
            )
          ),
          use.names = TRUE, fill = TRUE
        )
      }
      
      allACCs <- allACCs$Protein
      gnDT <- data.table()
      div <- length(allACCs) %% n.call # maximum number of values to call is 100
      mod <- length(allACCs) %/% n.call
      
      message(
        "Corresponding gene names of ", length(allACCs), " proteins will be called from uniprot in ",
        mod + 1, " iterations of maximum ", n.call, " calls."
      )
      
      if (mod != 0) {
        for (ix in 1:mod) {
          print(ix)
          protlist <- allACCs[seq((ix - 1) * (n.call) + 1, ix * n.call)]
          gnDT <- .calluniprot(protlist, gnDT)
        }
        if (div != 0) {
          protlist <- allACCs[seq((mod * n.call) + 1, length(allACCs))]
          gnDT <- .calluniprot(protlist, gnDT)
        }
      } else {
        protlist <- allACCs
        gnDT <- .calluniprot(allACCs, gnDT)
      }
      
      setnames(gnDT, uniprot.KEY, uniprot.input)
      if (toupper(uniprot.output) == "GENES") {
        gnDT <- gnDT[, Genes := tstrsplit(GENES, " ", fixed = TRUE)[1]][, -"GENES"] # first gene as the gene name
      }
    }
    
  }






  protDT. <- .lookupSplit(dpg = protDT., ref = gnDT, output = uniprot.output, input = uniprot.input)
  protDT.[, (uniprot.output) := do.call(paste, c(.SD, sep = split.char)), .SDcols = cols]
  protDT.[, (uniprot.output) := lapply(.SD, function(x) {
    gsub(paste0(split.char, "NA"), "", x)
  }), .SDcols = uniprot.output]
  GnPrtList <- protDT[, (cols) := NULL][, (uniprot.output) := protDT.[, get(uniprot.output)]][, -c("cntProt")]


  ## shorten gene names (if multiple genes are combined)
  ## example: Atp2b1;Atp2b2 --> ATP2b/1/2
  if (shorten.gene.names) {
    GnPrtList[, Genes.short := Genes]
    GnPrtList[Genes.short %like% split.char, Genes.short := sapply(Genes.short, function(x) .cutGeneName(x, split.char = split.char))]
    message("\n\nTrimmed gene names will be used in visualization!\nExample: 'Atp2b/1/2' instead of 'Atp2b1;Atp2b2'")
  }

  # if after shorting the gene name, there are genes with // or /// at the end!
  # (this happens when proteins in a master protein return the same gene! )
  dummy <- GnPrtList[Genes.short %like% "//", list(Protein, Genes.short)] %>% unique()
  if (nrow(dummy) > 0) {
    for (ix in c("////", "///", "//")) {
      dummy$gn <- gsub(ix, "", dummy$Genes.short)
    }
    GnPrtList <- merge(GnPrtList, dummy[, list(Protein, gn)], by = "Protein", all = TRUE)
    GnPrtList[!is.na(gn), Genes.short := gn][, gn := NULL]
  }

  # if different proteins map to same genes
  dup.gn <- GnPrtList[, uniqueN(Protein), by = Genes.short][V1 != 1, Genes.short]
  GnPrtList[Genes.short %in% dup.gn, Genes.short := paste0(Genes, "_Prt:", Protein)]

  MSqb2::char2fact(GnPrtList)
  MSqb2::char2fact(dat)


  if (save.annotation) {
    write.xlsx(GnPrtList,
      file = file.path(save.annotation.path, "AnnotationsTable.xlsx"),
      overwrite = TRUE
    )
  }

  dat <- merge(dat, GnPrtList[, .(Protein, Genes.short)], by = uniprot.input, all.x = TRUE)
  setnames(dat, "Genes.short", "Genes")


  dat[is.na(Genes) | Genes == "NA", Genes := factor(paste0("Gene: NA, Prt: ", get(uniprot.input)))]

  return(dat)
}



###################### merge multi column data with lookup table
.lookupSplit <- function(dpg, ref, input, output) {
  setkeyv(ref, cols = input)

  dpg[, rn := seq(nrow(dpg))] %>%
    # wide to long
    melt(., id.vars = "rn", variable.name = "inputGroup", value.name = "input") %>%
    .[!is.na(input)] %>%
    # Set key for the join
    setkey(., input) %>%
    ref[.] %>%
    # wide to long
    dcast(., rn ~ inputGroup, value.var = output) %>%
    .[, rn := NULL] %>%
    return(.)
}


###################### shorten long gene-family names!
##  Rab10;Rab1A;Rab3a;Rab14;Rab37;Rab1b;Rab4b;Rab39b;Rab8b;Rab43 ->
##  Rab /10/1A/3a/14/37/1b/4b/39b/8b/43
.cutGeneName <- function(xc, split.char = ";") {
  strarr <- strsplit(as.character(xc), split = split.char)[[1]]
  strarr %>% strsplit(., "") -> charlist

  if (length(strarr) > 1) {
    Reduce(intersect, charlist) -> chars
    if (length(chars) > 0) {
      charpos <- lapply(charlist, function(x) which(x %in% chars)) %>% Reduce(intersect, .)

      if (length(charpos) > 0) {
        ix <- which(diff(charpos, lag = 1) != 1)
        if (length(ix) > 0) charpos <- seq(ix)

        if (min(charpos) == 1) {
          charlist[[1]][charpos] %>% paste(., collapse = "") -> bs
          strarr <- paste0(bs, paste(gsub(bs, "/", strarr), collapse = ""))
        } else {
          strarr <- xc
        }
      } else {
        strarr <- xc
      }
    } else {
      strarr <- xc
    }
  } else {
    strarr <- xc
  }

  return(strarr)
}
