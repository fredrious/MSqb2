
msqb_config <- function(build.para,
                        config.file = NULL,
                        interactive = TRUE,
                        #### --------------- Filters
                        study.variable = "protein",
                        filter.by.quaninfo = "auto",
                        isolation.interference.cutoff = 75,
                        multi.features.method = "ionscore",
                        filter.singleshot.proteins = "ByPeptide",
                        min.intensity = 0.01,
                        reference.channel = NULL,
                        #### --------------- Normalization
                        normalization.method = "vsn",
                        normalize.bySubset = NULL,
                        #### --------------- Summarization
                        summarization.method = "tmp",
                        medianpolish.method = "column.eff",
                        medianpolish.byPool = FALSE,
                        #### ---------------Ftaction combination
                        fraction.collapse.method = "maxset",
                        #### --------------- remove batch effect
                        batch.corr.method = "none",
                        #### --------------- impute missing data
                        na.imputation.method = "none",
                        #### --------------- Filter metadata
                        ## (this will consequently filter the dataset)
                        filter.metadata = NULL,
                        #### --------------- Annotation
                        feature.annotation.source = "psm.file",
                        study.organism = NULL,
                        #### --------------- defining contrasts
                        pairwise.contrasts = TRUE,
                        manual.contrasts = NULL,
                        choose.contrasts = FALSE,
                        pairwise.denominator = NULL,
                        #### ---------------input data format (ProteomDiscoverer)
                        protein.column = "Master Protein Accessions",
                        peptide.seq.column = "Sequence",
                        peptide.annot.seq.column = "Annotated Sequence",
                        charge.column = "Charge",
                        FileID.column = "File ID",
                        quaninfo.column = "Quan Info",
                        isolation.interference.column = "Isolation Interference [%]",
                        ion.score.column = "Ions Score",
                        filename.column = "Spectrum File",
                        abundance.columns = "Abundance",
                        desc.column = "Master Protein Descriptions",
                        contaminant.column = "Contaminant",
                        marked.as.column = "Marked as",
                        #### --------------- control validating parameters/columns names
                        subsymDT = data.table(
                          symout = c("\\~", "\\+", "\\:", "\\|", "\\-", "\\(", "\\)", "\\/", " "),
                          symin = c("_", ".p", ".", ".", "_", ".", ".", ".", ".")
                        ),
                        ...) {

   ## assert main input para
   MSqb2:::.loggAssert(assertList(build.para))
   MSqb2:::.loggAssert(assertCharacter(config.file, len = 1, null.ok = TRUE))
   MSqb2:::.loggAssert(assertLogical(interactive, len = 1, null.ok = FALSE))


   ## import args from build parameters
   list2env(
      mget(
         c("Scripts.path", "Tables.path", "analysis.name", "Data.path",
           "add.date.tag", "interactive"),
         envir = as.environment(build.para)),
      envir = environment())


   ## list all Args manually passed to msqb_config function
   conf.args <- as.list(match.call(expand.dots = TRUE))
   conf.args[which(names(conf.args) == "")] <- NULL
   conf.args[which(names(conf.args) %in% c("config.file", "build.para"))] <- NULL
   conf.args[(conf.args == "F")] <- FALSE
   conf.args[(conf.args == "T")] <- TRUE


   ## importing config parameters from config file
   ## If no config.file provided, create one using the internal template
   conf.path <- .importConfigfile(
      conf.fl = config.file,
      Scripts.path = Scripts.path,
      analysis.name = analysis.name,
      add.date.tag = add.date.tag,
      whichConf = "workflow",
      conf.args = conf.args)


  ## after importing config parameters from config file, force all
  ## parameters in ellipsis
  .expandList(conf.args)


  ## add parameters which are going to be [possibly] modified to a second list
  conf.args2 <- list()



  #==================
  ## check if study must be done in Protein-level or in Peptide-level
  if (!isTRUE(check_choice(tolower(study.variable), c("protein", "peptide")))) {
    MSqb2:::.logg(FATAL, glue("'study.variable' must either be 'Protein' (for protein-leve analysis) ",
          "or 'Peptide' (for peptide-level analysis)."))
  }





  #==================
  MSqb2:::.loggAssert(assertCharacter(filter.by.quaninfo))
  MSqb2:::.loggAssert(assertNumber(isolation.interference.cutoff, null.ok = TRUE, lower = 0, upper = 100))
  MSqb2:::.loggAssert(assertNumber(min.intensity, null.ok = TRUE, lower = 0))
  MSqb2:::.loggAssert(assertCharacter(reference.channel, null.ok = TRUE))


  if (!isTRUE(check_choice(tolower(multi.features.method), c('ionscore', 'mean', 'median', 'max', 'all')))) {
     MSqb2:::.logg(FATAL, glue("'multi.features.method' must be one of the followings:\n",
                       "{glue_collapse(c('IonScore', 'mean', 'median', 'max', 'all'), sep = ', ')} \n"))
  }

  if (!isTRUE(check_choice(tolower(filter.singleshot.proteins), c('bypeptide', 'byfeature')))) {
     MSqb2:::.logg(FATAL, glue("'filter.singleshot.proteins' must be one of the followings:\n",
                       "{glue_collapse(c('bypeptide', 'byfeature'), sep = ', ')} \n"))
  }







  #==================
  ## check if Fraction combination method is among the available choices ----
  fraction.methods <- c("MAXset", "max", "median", "sum", "single", "all", "multi")
  if (!exists("fraction.collapse.method") ||
      !toupper(fraction.collapse.method) %in% toupper(fraction.methods)) {
      MSqb2:::.logg(WARN, glue(
        "Fraction Combination method must be one of the followings: \n",
        "{glue_collapse(fraction.methods, sep = ', ')} \n",
        "By default, 'fraction.collapse.method' will be set to 'MAXset'. \n",
        "Alternatively you can edit the config file and rerun msqb_config function."
      ))
    conf.args2[["fraction.collapse.method"]] <- "MAXset"
  }





  #==================
  ## check if normalization method is among the available choices ----
  normalization.methods <- c("vsn", "median", "quantile")
  if (!exists("normalization.method") ||
      !toupper(normalization.method) %in% toupper(normalization.methods)) {
       MSqb2:::.logg(WARN, glue(
        "Normalization method must be one of the followings: \n",
        "{glue_collapse(normalization.methods, sep = ', ')} \n",
        "By default, 'normalization.method' will be set to 'vsn'. \n",
        "Alternatively you can edit the config file and rerun msqb_config function."
      ))
    conf.args2[["normalization.method"]] <- "vsn"
  }

  if (!isTRUE(check_character(normalize.bySubset, null.ok = TRUE, len = 3))) {
     MSqb2:::.logg(ERROR, glue("'normalize.bySubset' must be a list with three character elements:\n",
                       "Element 1 must be the name of a column in the dataset that is to be used ",
                       "for subsetting.\n",
                       "Element 2 is the character string in the selected column that is to be ",
                       "used for row indexing.\n",
                       "Element 3 must either be 'exact', meaning that matching must be exact, or ",
                       "'partial', which allows partial matching.\n"))
  }

  if (!is.null(normalize.bySubset) &
      !isTRUE(check_choice(normalize.bySubset[3],choices = c("exact", "partial") ))) {
     MSqb2:::.logg(ERROR, glue("In 'normalize.bySubset', element 3 must either be 'exact', ",
     "meaning that matching must be exact, or 'partial', which allows partial matching.\n"))
  }







  #==================
  ## check parameters passed to filter metadata ----
  if (!is.null(filter.metadata)) {
    filter.metadata <- eval(filter.metadata)
    conf.args2[["filter.metadata"]] <- filter.metadata
  }






  #==================
  ## check parameters for missing value imputation ----
  imp.methods <- c("knn", "QRILC", "MinDet", "MinProb", "Zero", "mle", "Hybrid", "none")
    if (!exists("na.imputation.method", mode = "character") ||
      !toupper(na.imputation.method) %in% toupper(imp.methods)) {
         MSqb2:::.logg(WARN, glue(
          "Missing value imputation method must be one of the followings: \n",
          "(From 'imputeLCMD' package:) \n",
          "{glue_collapse(imp.methods, sep = ', ')} \n",
          "By default, 'na.imputation.method' will be set to 'none'."
        ))
      conf.args2[["na.imputation.method"]] <- "none"
    }





  #==================
  ## check advanced mode ----
  if (!exists("advance.mode", mode = "logical" )) {
    conf.args2[["advance.mode"]] <- FALSE
  } else if (advance.mode & !exists("model.formula", mode = "character")) {
    MSqb2:::.logg(ERROR, glue(
       "advance.mode is TRUE.\n",
       "model.formula must be provided! example: '~ 0 + Condition'"))
    model.formula <- as.character(deparse(substitute(model.formula)))
  }





  #==================
  ## check model.formula (to be passed to limma) ----
  if (exists("model.formula") &
      !exists("advance.mode", mode = "logical"  )) {
     MSqb2:::.logg(WARN, glue(
        "model.formula is provided. Automatic experimental-design detection will be deactivated.",
        "Analysis mode: advance"))
        conf.args2[["advance.mode"]] <- TRUE
    if (!is.character(model.formula) ) { #if expression
          model.formula <- as.character(deparse(substitute(model.formula)))
    }
  }





  #==================
  ## check for batch effect parameters ----
  rmbtch.methods <- c("limma", "svn", "none")
  if (!exists("batch.corr.method") ||
      !toupper(batch.corr.method) %in% toupper(rmbtch.methods)) {
       MSqb2:::.logg(WARN, glue(
        "Batch effect removal method must be one of the followings: \n",
        "{glue_collapse(rmbtch.methods, sep = ', ')} \n",
        "By default, 'batch.corr.method' will be set to 'none'."
      ))
    conf.args2[["batch.corr.method"]] <- "none"
  }


  if (toupper(batch.corr.method) == "LIMMA") {
     MSqb2:::.logg(TRACE, glue(
        "\n!!! When batch correction method is set to 'LIMMA', the batch-corrected data will ",
        "only be used for data visualization, like in PCA and heatmaps. As recommended in the description of ",
        "the function 'removeBatchEffect' from LIMMA package, the batch factors will be included ",
        "in the linear model, if not already included by the user!",
        "\nExample: if model formula is '~ 0 + Condition' and the batch factor is 'Pool' ",
        "the model formula will be set to '~ 0 + Condition + Pool'.\n",
        "Two batch parameters can be used for batch correction. By default, the first ",
        "two factors after Condition in the design formula will be used."
      ))
  } 

  # if (toupper(batch.corr.method) == "COMBAT") {
  #      MSqb2:::.logg(WARN, glue(
  #       "For batch correction using method 'ComBat', no missing values in the data is ",
  #       "permitted. Missing value imputation (if not already set by user) will automatically ",
  #       "be performed in advanced, using imputation method 'KNN'!\n",
  #       "To skip batch correction, set batch.corr.method to 'none'."))
  # 
  #     if (!exists("na.imputation.method") && !"na.imputation.method" %in% names(conf.args2) ){
  #       conf.args2[["na.imputation.method"]] <- "KNN"
  #     }
  # }







  #==================
  ## check parameters for peptide to protein summarization ----
  pep2prot.methods <- c("TMP", "median", "mean", "sum")

  if (!exists("summarization.method", mode = "character") ||
      !toupper(summarization.method) %in% toupper(pep2prot.methods)) {
       MSqb2:::.logg(WARN, glue(
        "peptide to protein summarization method must be one of the followings: \n",
        "{glue_collapse(pep2prot.methods, sep = ', ')} \n",
        "By default, 'summarization.method' will be set to 'TMP' (Tukey Median Polish)."
      ))
    conf.args2[["summarization.method"]] <- "TMP"
  }

  if (toupper(summarization.method) == "TMP") {
    if (!exists("medianpolish.method") ||
      !toupper(medianpolish.method) %in% toupper(c("column.eff", "residual"))) {
       MSqb2:::.logg(WARN, glue(
          "When peptide to protein summarization method is set to 'TMP' (Tukey Median Polish), ",
          "parameter medianpolish.method must be set to either of 'column.eff' and 'residual'. ",
          "By default, 'medianpolish.method' will be set to 'column.eff'."
       ))
      conf.args2[["medianpolish.method"]] <- "column.eff"
    }

    if (!exists("medianpolish.byPool", mode = "logical") ||
        !is.logical(medianpolish.byPool)) {
      conf.args2[["medianpolish.byPool"]] <- FALSE
    }
  }






  #==================
  ## check parameters for gene annotation  ----
  # check where to look for gene annotations
  # annt.src <- c("annotationdbi", "uniprot", "none")
  annt.src <- c("psm.file", "annotationdbi", "uniprot.web", "uniprot.file", "none")
  do.annotation <- TRUE # default

  if (isTRUE(check_choice(tolower(feature.annotation.source), annt.src))) {
  
    if (tolower(feature.annotation.source) == "none") {
      do.annotation <- FALSE
    } else if (tolower(feature.annotation.source) %in% c("annotationdbi", "uniprot")) {
      
      orgs <- c("human", "mouse", "yeast")
      if (!exists("study.organism")) {
        MSqb2:::.logg(ERROR, glue("'feature.annotation.source' is set to {feature.annotation.source}.\n",
                                 "'study.organism' must also be defined and set to either of the followings:\n",
                                 "{glue_collapse(orgs, sep = ', ')}"))
      }
    } else if (tolower(feature.annotation.source) == "psm.file") {
      feature.annotation.source <- build.para$measurements.file
    }
  } else {
    conf.args2[["feature.annotation.source"]] <-  
      list("feature.annotation.source" = feature.annotation.source) %>% 
      .check.file(., dpth = Data.path, toData = interactive, interactive = interactive) %>% 
      unlist() %>% unname()
  }
  
  





  #==================
  ## check parameters for statistical analysis with Limma  ----

  if (!exists("limma.blocking.parameter", mode = "character") ||
    !exists("limma.blocking.parameter", mode = "NULL")) {
    conf.args2[["limma.blocking.parameter"]] <- NULL
  }

  # limma.blocking.parameter
  # if (exists("limma.blocking.parameter", mode = "character")) {
  #   techrep <- limma.blocking.parameter
  #   conf.args2[["techrep"]] <- limma.blocking.parameter
  # }
  #
  # # limma.trend
  # if (!exists("limma.trend", mode = "logical")) {
  #   MSqb2:::.logg(WARN, glue(
  #     "Parameter 'limma.trend' for limma statistical test ",
  #     "must be a logical value. By default it will be set to FALSE."
  #   ))
  #   conf.args2[["limma.trend"]] <- TRUE
  # }

  # pairwise.contrasts
  if (!exists("pairwise.contrasts", mode = "logical")) {
    MSqb2:::.logg(WARN, glue(
      "Parameter 'pairwise.contrasts' for limma statistical test ",
      "must be a logical value. By default it will be set to TRUE."
    ))
    conf.args2[["pairwise.contrasts"]] <- TRUE
  }

  # manual contrast
  if (exists("manual.contrasts") & storage.mode(manual.contrasts) == "language") {
    manual.contrasts <- eval(manual.contrasts)
    conf.args2[["manual.contrasts"]] <- manual.contrasts
  }
  if (!exists("manual.contrasts", mode = "character") &&
      !exists("manual.contrasts", mode = "NULL")) {
    MSqb2:::.logg(WARN, glue(
      "Parameter 'manual.contrasts' (manually defined pairwise contrasts) must either be NULL or a ",
      "'character' vector. By default it will be set to NULL."
    ))
    conf.args2[["manual.contrasts"]] <- NULL
  }











  ## remove items of conf.args already exist in conf.args2
  conf.args[names(conf.args) %in% names(conf.args2)] <- NULL
  ## now merge the two conf.args
  allargs <- c(conf.args2, conf.args)
  allargs$build.para <- NULL



  ## after parameter check, save changes to config parameters in config.file
  if (length(allargs) > 0) {
    for (i in 1:length(allargs)) {
      assign(names(allargs)[i], allargs[[i]])

      txtin <- paste(names(allargs)[i], "=", paste(deparse(allargs[[i]]), collapse = ""))
      code <- readLines(conf.path) %>% gsub("<-", "=", .)
      nl <- which(lapply(strsplit(gsub(" ", "", code), split = "=", ), "[", 1) == names(allargs)[[i]])
      if (length(nl) > 1) {
        MSqb2:::.logg(FATAL, glue(
          "Multiple matches for the followings parameters: \n",
          "{glue_collapse(names(allargs)[[i]], sep = ', ')} \n",
          "Edit the config file and rerun the msqb_config."
        ))
      } else if (length(nl) == 1) {
        txtout <- code[nl]
        code.edited <- gsub(txtout, txtin, code, fixed = TRUE)
        writeLines(code.edited, conf.path)
        rm(code.edited)
      } else if (length(nl) == 0) {
        write(txtin, file = conf.path, append = TRUE)
      }
    }
  }

  if (interactive) {

    MSqb2:::.logg(TRACE, glue(
      "Do you want to view/edit the config file? \n",
      "(If yes, save changes after editting and before ",
      "starting the workflow."))

   res <- select.list(c("yes", "no"))
   if (res == "yes") file.edit(conf.path)
   Sys.sleep(0.1)
  }

  ## save configuration parameters in the Analyses/Tables directory.
  ## The file will be tagged by date and time for reproducibility purposes.
  if (!dir.exists(file.path(Tables.path, "Config_Parameters"))) {
    dir.create(file.path(Tables.path, "Config_Parameters"))
  }

  paste0("_", format(Sys.time(), "%Y%m%d")) %>%
    paste0("ConfigParameters_", ., ".txt") %>%
    file.path(file.path(Tables.path, "Config_Parameters"), .) %>%
    file.copy(from = conf.path, to = ., overwrite = FALSE)


  suppressWarnings(
     rm(txtin, txtout, nl, code, i, res, fi, fls,
       normalization.methods, fraction.methods, allargs, conf.args, conf.args2,
       imp.methods, Scripts.path, analysis.name, add.date.tag, gene.annotation.file
  ))

  conf.vars <- mget(ls(name = environment()), envir = environment()) %>%
    lapply(., eval)

  return(conf.vars)
}
