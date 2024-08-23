#' MSqb Workflow for Proteomics Data Analysis
#'
#' This function orchestrates the entire workflow for MS proteomics data analysis, including data reformatting, normalization, batch effect removal, protein summarization, differential expression analysis, and more.
#'
#' @param build.para A list of parameters required to build the workflow environment (see `MSqb_build()` function).
#' @param config.para A list of configuration parameters for the workflow. If NULL, these will be read from the `config.file` (see `MSqb_config()` function).
#' @param config.file The path to the configuration file. This file will be used if `config.para` is not provided.
#' @param ... Additional parameters passed to internal functions.
#'
#' @return An environment containing the results of the workflow, including processed data tables and statistical results.
#' @export
#' @importFrom glue glue
#'
#' @examples
#' \dontrun{
#' result_env <- msqb_workflow(
#'   build.para = list(param1 = value1, param2 = value2),
#'   config.file = "config.R"
#' )
#' }
msqb_workflow <- function(build.para,
                          config.para = NULL,
                          config.file = NULL,
                          ...) {



   #==================
   #### make new env. to store results and parameters
   res.env <- new.env()
   list2env(list("config.para" = config.para, "build.para" = build.para), envir = res.env)




   #==================
   #### i) read config parameters ----
   ## import config parameters (either from config.file, config.para or ellipsis!)
   ## and build parameters to msqb_workflow environment.
   readConfigPara(
      build.para = build.para,
      config.para = config.para,
      config.file = config.file,
      config.type = "wf")



   #==================
   #### ii) read metadata ----
   metadata0 <- read.file(file = metadata.file, sheet = metadata.file.sheet) %>% unique()
   ## Add BioRep if not already in the metadata
   exp.design <- check_design(metadata0)
   ## check if parameters in the model formula match column names of the metadata
   model.formula <- .match_ModelFormula_metadata(metadata0, model.formula)




   #==================
   #### ii) read measurement data ----
   msDT <- data.table()
   metadata <- data.table()
   annotdata <- data.table()

   ## File.Id MUST be different in multiple PSM files if input is from PD
   for (id in seq_along(measurements.file)) {

    msDT.dummy <- read.file(file = measurements.file[id],
                             sheet = measurements.file.sheet)

    if (all(c(filename.column, FileID.column) %in% names(msDT.dummy))) {
       msDT.dummy[, .(filename.column = uniqueN(get(filename.column))), by = FileID.column] %>%
          { .$filename.column > 1 } %>%
          any(.) %>%
          { if(.) { .logg(ERROR, glue(
             "Some {FileID.column}s match by multiple raw files! ",
             "In case of multiple input files, make sure that the {FileID.column}s are unique!"))
             }
          }
    }



    #==================
    ## 1. Reformat data ----
    ## !output is a list!
    ## Two objects will be exported to the Global env.: ms.software and MSmethod
    msDT.dummy <- reformat_PD_data(
      feat.dt = msDT.dummy,
      metadata = metadata0,
      protein.column = protein.column,
      peptide.seq.column = peptide.seq.column,
      peptide.annot.seq.column = peptide.annot.seq.column,
      charge.column = charge.column,
      FileID.column = FileID.column,
      quaninfo.column = quaninfo.column,
      isolation.interference.column = isolation.interference.column,
      abundance.columns = abundance.columns,
      ms.software = ms.software,
      ion.score.column = ion.score.column,
      desc.column = desc.column,
      contaminant.column = contaminant.column,
      marked.as.column = marked.as.column,
      subsymDT = subsymDT
    )

    metadata <- rbind(metadata, msDT.dummy$metadata) %>% unique()# categorical vars.
    annotdata <- rbind(annotdata, msDT.dummy$annotdata) # protein/gene annotation
    msDT <- rbind(msDT, msDT.dummy$feat.dt) # main DT
    rm(msDT.dummy)
  }
   # store in res environemnt
   list2env(list("annotdata" = annotdata), envir = res.env)




  ## .........................
  ## 2. Remove Sample(s) ----
  ## if any sample is to be removed! (a channel, replicate or fraction)
  ## This is usually done after visual inspection of the plots after the first run.
  if (!is.null(filter.metadata)) {
    metadata <- filter_metadata(metadata, filters = eval(filter.metadata))
    msDT <- merge(msDT, metadata, by = names(metadata), all.y = TRUE)
  }
  # store in res environemnt
  list2env(list("metadata" = metadata), envir = res.env)





  ## ....................
  ## 3. Filter data ----
  ## Two objects will be added to res.env: BasicStats and Isol.Interference
  msDT <- clean_feature_data(
    dat = msDT,
    metadata = metadata,
    filter.by.quaninfo = filter.by.quaninfo,
    isolation.interference.cutoff = isolation.interference.cutoff,
    collapse_psm_method = collapse_psm_method,
    min.intensity = min.intensity
  )


  # store in res environemnt
  list2env(list(
    "FeatureDT.cleaned" = msDT$msDT,
    "BasicStats" = msDT$BasicStats,
    "Isol.Interference" = msDT$Isol.Interference
  ), envir = res.env)
  msDT <- msDT$msDT






  ## .........................................
  ## 4. Make summary tables for qc plots ----
  QCtables <- mkQCtables(
    msDT = msDT,
    metadata = metadata,
    ms.software = ms.software,
    Tables.path = Tables.path
  )
  ## add QC tables to res.env
  list2env(QCtables, envir = res.env)





#list2env(list("PeptideDT.normalised" = msDT), envir = res.env)
  ## .............................
  ## 5. Fraction combination ----
  if ("Fraction" %in% names(msDT)) {
    msDT <- collapse_fractions(
      dat = msDT,
      group = "Pool",
      method = fraction.collapse.method,
      path = Tables.path
    )
    metadata[, Fraction := NULL]
    list2env(list("PrtPepCount_fraction" = msDT$PrtPepCount_fraction), envir = res.env)
    msDT <- msDT$dt.fcomb
    metadata <- metadata[, Filename := do.call(paste, c(.SD, sep = "_")), .SDcols = intersect(names(msDT), c("Pool", "Channel"))] %>% unique(.) %>% char2fact()
    msDT <- msDT[, Filename := do.call(paste, c(.SD, sep = "_")), .SDcols = intersect(names(msDT), c("Pool", "Channel"))] %>% unique(.)
  }





  ## ......................
  ## 6. Normalization ----
  msDT <- normalise_data(
    dat = msDT,
    method = normalization.method,
    val = "Intensity",
    calib = "affine", # if method vsn
    normalize.bySubset = normalize.bySubset,
    path = Tables.path
  )
  # store in res environemnt
  list2env(list("PeptideDT.normalised" = msDT), envir = res.env)






  ## ............................
  ## 7. impute missing data ----
  if (na.imputation.method != "none" | toupper(batch.corr.method) == "COMBAT") {
    msDT <- impute_missing_values(
      msDT,
      val = "Abundance",
      col.p = "Filename",
      row.p = intersect(c("Protein", "Peptide", "Feature"), names(msDT)),
      na.imputation.method = na.imputation.method,
      hybrid.mar = "KNN", # if method "Hybrid"
      hybrid.mnar = "QRILC", # if method "Hybrid"
      n_knn = n_knn, # knn imputation
      tune.sigma = 0.5, # QRILC, MinProb imputation
      q = 0.01
    )

    if (length(setdiff(names(metadata), names(msDT))) > 0) {
       msDT <- merge(msDT, metadata,
                     by = intersect(names(msDT), names(metadata)),
                     all.y = TRUE)
    }
  }





  ## ............................
  ## 8. remove batch effect ----
  ## currently only limma
  if (batch.corr.method != "none") {
    PeptideDT.batchCorrected <- remove_batch_effect(
      dat = msDT,
      metadt = metadata,
      batch.corr.method = batch.corr.method,
      val = "Abundance",
      model.formula = model.formula
    )

    if (toupper(batch.corr.method) == "LIMMA") {
      # no change in input data, just for visualisation
      list2env(list("PeptideDT.batchCorrected" = PeptideDT.batchCorrected), envir = res.env)
    } else {
      msDT <- PeptideDT.batchCorrected
      rm(PeptideDT.batchCorrected)
    }
  }





  ## ......................
  ## 9. summarization ----
  msDT <- protein_aggregation(
    dat = msDT,
    summarization.method = summarization.method,
    medianpolish.method = medianpolish.method,
    medianpolish.byPool = medianpolish.byPool,
    block = "Pool",
    var = "Abundance",
    groupby = capitalize_first_char(study.variable) ,
    Tables.path = Tables.path
  )

  msDT <- merge(msDT, metadata, by = intersect(names(msDT), names(metadata)), all.x = TRUE)
  if (exists("annotdata") &&!is.null(annotdata) && nrow(annotdata) > 1) {
    msDT <- merge(msDT, annotdata, by = intersect(names(msDT), names(annotdata)), all.x = TRUE)
  }
  # store in res environment
  list2env(list("ProteinDT" = msDT), envir = res.env)



  ## store PSM count per protein, to be passed to limma trend
  msDT <- add_psm_count(msDT, param = capitalize_first_char(study.variable))
  msDT <- merge(msDT, metadata, by = intersect(names(msDT),names(metadata)))


  ## .........................
  ## 10. Gene annotation ----
  # if (do.annotation) {
  #
  #   msDT <- annotation_query(
  #     dat = msDT,
  #     feature.annotation.source = feature.annotation.source, # "uniprot", "uniprot.file", "annotation.file"
  #
  #     # if source "uniprot":
  #     uniprot.input = "Protein",
  #     uniprot.output = "Genes",
  #     uniprot.KEY = "UNIPROTKB",
  #     uniprot.species.name = uniprot.species.name,
  #     uniprot.taxId = uniprot.taxId,
  #     n.call = 99,
  #
  #     # if source "annotation.file":
  #     # annotation.file = "221128_P-481-CL_Pool1_2_3 221202 prot.xlsx",
  #     shorten.gene.names = TRUE, # shorten the gene and protein names?
  #     split.char = ";",
  #     save.annotation = TRUE, # only if annot.sourc == "uniprot"
  #     save.annotation.path = Tables.path,
  #     Data.path = Data.path
  #   )
  #   # store in res environemnt
  #   list2env(list("ProteinDT" = msDT), envir = res.env)
  # }





  ## ...................
  ## 12. if batch corrected data exist ----
  ## repeat summarization and annotation
  if (exists("PeptideDT.batchCorrected")) {

    ProteinDT.batchCorrected <- protein_aggregation(
      dat = PeptideDT.batchCorrected,
      summarization.method = summarization.method,
      medianpolish.method = medianpolish.method,
      medianpolish.byPool = medianpolish.byPool,
      block = "Pool",
      var = "Abundance",
      groupby = capitalize_first_char(study.variable),
      Tables.path = Tables.path
    )
    ProteinDT.batchCorrected <- merge(ProteinDT.batchCorrected,
                                      unique(metadata[, !"Fraction"]),
                                      by = intersect(names(ProteinDT.batchCorrected),
                                                     names(metadata)), all.x = TRUE)
    if (!is.null(annotdata)) {
      ProteinDT.batchCorrected <- merge(ProteinDT.batchCorrected,
                                        annotdata,
                                        by = intersect(names(ProteinDT.batchCorrected),
                                                       names(annotdata)), all.x = TRUE)
    }

    # if (do.annotation) {
    #   char2fact(ProteinDT.batchCorrected)
    #   char2fact(msDT)
    #   ProteinDT.batchCorrected <- merge.data.table(ProteinDT.batchCorrected, msDT[, -"Abundance"])
    # }
    # store in res environemnt
    list2env(list("ProteinDT.batchCorrected" = ProteinDT.batchCorrected), envir = res.env)
  }




  ## ...................
  ## 11. fit Limma ----
  limma.fit <- DE_stat_func(
    dat = msDT,
    sample = "Filename",
    model.formula = model.formula, # model.formula,
    variable = capitalize_first_char(study.variable),
    value = "Abundance",
    techrep = "TechRep",
    pairwise.contrasts = pairwise.contrasts,
    pairwise.denominator = pairwise.denominator,
    manual.contrasts = manual.contrasts,
    choose.contrasts = choose.contrasts,
    make.complex.contrasts = make.complex.contrasts,

    # for repeated measurements and/or technical Reps
    blocking.parameter = limma.blocking.parameter,
    limma.trend = limma.trend, # check limma plotSA!
    limma.robust = limma.robust,
    limma.winsor.tail.p = FALSE,
    time.course = FALSE,
    removeTotalMissing = TRUE,
    reverse.contrasts = FALSE,
    rmSingleShotProteins = TRUE,
    add.vs.inContrast = FALSE,
    path = Tables.path
  )
  # store in res environemnt
  list2env(list(
    "fit_residuals" = limma.fit$fit_residuals,
    "fit_eb" = limma.fit$fit_eb,
    "fit0" = limma.fit$fit0
  ),
  envir = res.env
  )
  limma.fit <- limma.fit$fit_eb






  ## ...........................
  ## 12. generate topTable ----
  statLists <- DE_stat_table(
    fit_eb = limma.fit,
    dat = msDT,
    variable = "Protein",
    path = Tables.path,
    adjust = "BH",
    Gene.name = ifelse("Genes" %in% names(msDT), "Genes", "Protein"),
    save = TRUE,
    file.name = "limma.stat_output"
  )
  # store in res environemnt
  list2env(list("topList" = statLists[order(Comparison, rank.id)]), envir = res.env)




  ## save session info and image ----
  writeLines(
    text = capture.output(sessionInfo()),
    con = file.path(Analysis.path, paste0(prefix, "_SessionInfo_", suffix, ".txt"))
  )
  save(res.env, file = file.path(Analysis.path, paste0(prefix, "_Results_Image_", suffix, add.date.tag, ".Rdata")))

  return(res.env)
}


