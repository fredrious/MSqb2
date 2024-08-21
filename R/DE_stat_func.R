#' Perform Differential Expression Analysis Using LIMMA
#'
#' This function performs differential expression analysis on MS data using limma (Linear Models for Microarray Data) package.
#' It supports a wide range of experimental designs, including handling of technical replicates, blocking factors, and the application of contrasts
#' to identify differentially expressed proteins or peptides.
#'
#' @param dat A `data.table` containing the MS data. The table should include columns for samples, variables (e.g., proteins),
#'   and values (e.g., abundances).
#' @param sample A character string specifying the column name that contains sample identifiers. Default is `"SampleID"`.
#' @param variable A character string specifying the column name that contains the variables of interest (e.g., `"Protein"`). Default is `"Protein"`.
#' @param value A character string specifying the column name that contains the values to be analyzed (e.g., `"Abundance"`). Default is `"Abundance"`.
#' @param model.formula A character string specifying the model formula to be used in the analysis. Default is `"~ 0 + Condition"`.
#' @param sva.obj An optional `sva` object for batch correction. Default is `NULL`.
#' @param techrep A character string specifying the column name that contains technical replicates. If provided, these will be used as blocking factors. Default is `NULL`.
#' @param limma.trend Logical or character. If `TRUE`, trend fitting is applied; if a character, it specifies the column name to be used for trend fitting. Default is `TRUE`.
#' @param limma.robust Logical. If `TRUE`, robust fitting is applied in the eBayes step. Default is `TRUE`.
#' @param limma.winsor.tail.p Logical or numeric. If `FALSE`, no winsorization is applied; if numeric, specifies the tail probability for winsorization. Default is `FALSE`.
#' @param time.course Logical. If `TRUE`, a time course analysis is performed using spline regression. Default is `FALSE`.
#' @param spline.df Integer. Degrees of freedom for spline regression in time course analysis. Default is `4`.
#' @param removeTotalMissing Logical. If `TRUE`, rows with all missing values are removed from the analysis. Default is `TRUE`.
#' @param pairwise.contrasts Logical. If `TRUE`, all pairwise contrasts between conditions are computed. Default is `TRUE`.
#' @param pairwise.denominator Character vector. Specifies the condition levels to be used as denominators in pairwise contrasts. Default includes `"wt"`, `"wildtype"`, `"control"`, etc.
#' @param manual.contrasts A list of manual contrasts to be applied. Default is `NULL`.
#' @param choose.contrasts Logical. If `TRUE`, allows manual selection of contrasts from a list. Default is `FALSE`.
#' @param make.complex.contrasts Logical. If `TRUE`, prompts the user to create complex contrasts. Default is `FALSE`.
#' @param blocking.parameter A character string specifying the column name for blocking factors, typically technical replicates. Default is `NULL`.
#' @param reverse.contrasts Logical. If `TRUE`, reverses the direction of the contrasts. Default is `FALSE`.
#' @param showWarnings Logical. If `TRUE`, warnings are shown during execution. Default is `FALSE`.
#' @param prefix A character string to prefix the output filenames. Default is `NULL`.
#' @param suffix A character string to suffix the output filenames. Default is `NULL`.
#' @param rmSingleShotProteins Logical. If `TRUE`, proteins with only one peptide are removed from the analysis. Default is `TRUE`.
#' @param add.vs.inContrast Logical. If `TRUE`, the contrasts are labeled with "vs." between conditions. Default is `TRUE`.
#' @param path A character string specifying the output directory for saving results. Default is the current working directory.
#' @param ... Additional arguments passed to the underlying LIMMA functions.
#'
#' @return A list containing:
#' \itemize{
#'   \item `"fit_eb"`: The fitted model after applying eBayes.
#'   \item `"fit_residuals"`: The residuals of the fit.
#'   \item `"fit0"`: The initial fitted model before applying contrasts.
#' }
#'
#' @details The `DE_stat_func` function is designed to handle complex experimental designs typical in mass spectrometry-based proteomics studies.
#' It supports:
#'
#' \itemize{
#'   \item **Technical Replicates**: When `techrep` is provided, the function accounts for repeated measures by using the `duplicateCorrelation` method from LIMMA.
#'   \item **Time Course Analysis**: If `time.course` is `TRUE`, the function performs spline regression to model the effect of time.
#'   \item **Contrast Definitions**: Users can define contrasts manually or choose from pairwise contrasts. Complex contrasts, such as average or delta-delta contrasts, can also be generated.
#' }
#'
#' This function is suitable for both differential expression analysis in standard proteomics experiments and more complex designs involving time courses or technical replicates.
#' When working with experiments that have potential batch effects, users can incorporate surrogate variable analysis (`sva.obj`) to correct for these.
#'
#' @importFrom limma lmFit makeContrasts eBayes duplicateCorrelation contrasts.fit
#' @importFrom splines ns
#' @importFrom stats model.matrix as.formula
#' @importFrom glue glue
#'
#' @export
DE_stat_func <- function(dat,
                      sample = "SampleID",
                      variable = "Protein",
                      value = "Abundance",
                      model.formula = "~ 0 + Condition",
                      sva.obj = NULL,
                      techrep = NULL,
                      limma.trend = TRUE,
                      limma.robust = TRUE,
                      limma.winsor.tail.p = FALSE,
                      time.course = FALSE,
                      spline.df = 4,
                      removeTotalMissing = TRUE,
                      pairwise.contrasts = TRUE,
                      pairwise.denominator = c("wt", "wildtype", "wild-type", "cntrl", "control", "ctl", "untreated"),
                      manual.contrasts = NULL,
                      choose.contrasts = FALSE,
                      make.complex.contrasts = FALSE,
                      blocking.parameter = NULL, # column corresponding to technical replicate
                      reverse.contrasts = FALSE,
                      showWarnings = FALSE,
                      prefix = NULL,
                      suffix = NULL,
                      rmSingleShotProteins = TRUE,
                      add.vs.inContrast = TRUE,
                      path = getwd(), ...) {
  # function body...
}

DE_stat_func <- function(data,
                     sample = "SampleID",
                     variable = "Protein",
                     value = "Abundance",
                     model.formula = "~ 0 + Condition",
                     sva.obj = NULL,
                     techrep = NULL,
                     limma.trend = TRUE,
                     limma.robust = TRUE,
                     limma.winsor.tail.p = FALSE,
                     time.course = FALSE,
                     spline.df = 4,
                     removeTotalMissing = TRUE,
                     pairwise.contrasts = TRUE,
                     pairwise.denominator = c("wt", "wildtype", "wild-type", "cntrl", "control", "ctl", "untreated"),
                     manual.contrasts = NULL,
                     choose.contrasts = FALSE,
                     make.complex.contrasts = FALSE,
                     blocking.parameter = NULL, # column corresponding to technical replicate
                     reverse.contrasts = FALSE,
                     showWarnings = FALSE,
                     prefix = NULL,
                     suffix = NULL,
                     rmSingleShotProteins = TRUE,
                     add.vs.inContrast = TRUE,
                     path = getwd(), ...) {

  # suppress warning
  oldw <- getOption("warn")
  if (!showWarnings) options(warn = -1)


  # if contrasts are not defines
  if (!pairwise.contrasts &
      !choose.contrasts &
      !make.complex.contrasts &
      !time.course &
      is.null(manual.contrasts)) {
    choose.contrasts <- TRUE
  }

  # contrasts will be stored in cnts
  cnts <- NULL
  repeated.measures <- FALSE


  if (rmSingleShotProteins & "Peptide.count" %in% names(data)) data <- data[Peptide.count != 1, ]

  ## check if data contains techrep.
  ## If it does, force it to be used as blocking factor
  if (!is.null(techrep)) {
    if (techrep %in% names(data)) {
      message(
        "Data contains technical replicates. ",
        "These will be passed to LIMMA as 'blocking' parameters!"
      )
      repeated.measures <- TRUE
      if (!is.null(blocking.parameter)) {
        if (blocking.parameter != techrep) {
          warning(
            "The 'blocking parameter' is different from 'technical replicate'! ",
            "The 'blocking.parameter' will be changed to 'techrep' by default!"
          )
          blocking.parameter <- techrep
        }
      } else {
        blocking.parameter <- techrep
      }
    }
  }



  if (is.null(blocking.parameter) & repeated.measures) {
    warning(
      "The parameter 'repeated.measures' is set to TRUE ",
      "but the blocking parameter is not provided!",
      "By default 'repeated.measures' will be set to FALSE!"
    )
    repeated.measures <- FALSE
  }

  if (!is.null(blocking.parameter) & !repeated.measures) {
    warning(
      "The blocking parameter is provided while the parameter 'repeated.measures' is FALSE!",
      "By default 'repeated.measures' will be set to TRUE!"
    )
    repeated.measures <- TRUE
  }


  # fit parameters
  fit.para <- c(all.vars(as.formula(model.formula), unique = TRUE), blocking.parameter)

  if (any(!fit.para %in% names(data))) {
    stop(
      "Some parameters used in design formula not match the ",
      "column names of input data!\n", model.formula
    )
  }

  # make reference (target) table
  reftb <- data[, c(..sample , ..fit.para)] %>% unique %>% setorderv(., fit.para)

  if (any(apply(reftb, 2, uniqueN) == 1)) {
    stop(
      "Parameter ", names(which(apply(reftb, 2, uniqueN) == 1)), " used in design formula has only 1 level!",
      "contrasts can be applied only to factors with at least 2 levels!"
    )
  }



  ############################################
  # time course analysis (spline regression limma)
  if (time.course) {
    if (!"time" %in% tolower(names(data))) {
      MSqb2:::.logg(FATAL, glue("'time.course' is set to TRUE but column 'time' does not exist in metadata. ",
                               "Fitting data considering many time points per group can not be performed!"))
    }


    if (pairwise.contrasts |  !is.null(manual.contrasts) | make.complex.contrasts) {
      MSqb2:::.logg(WARN, glue("No pairwise contrast will be generated in time course analysis."))
      pairwise.contrasts <- make.complex.contrasts <- add.vs.inContrast <- FALSE
      manual.contrasts <- NULL
    }
    cnts <- NULL # there won't be any contrast


    time.para <- grep("time", names(data), ignore.case = TRUE, value = TRUE)
    ctime <- data[, c(..sample , ..fit.para, ..time.para)] %>% unique %>% setorderv(., fit.para) %>%
      .[, ..time.para] %>% unlist(., use.names = FALSE) %>% as.factor() %>% as.numeric()

    # if time column is not numeric, try to convert it to numeric. If not possible, stop!
    # if (!is.numeric(ctime)) {
    #   MSqb2:::.logg(FATAL, glue("In 'time course' analysis, column {time.para} must be numeric. ",
    #                            "It is not and can not be converted to class numeric!"))
    # }


    # fit spline: look at Limma guide, section 9.6.2 --> Many time points
    library(splines)
    Bspline.mat <- splines::ns(ctime, df = spline.df)
    # update model formula
    model.formula %<>% gsub(fit.para[1], glue(fit.para[1], " * Bspline.mat"), .)
  }



  ############################################
  # define design matrix
  design.mat <- model.matrix(as.formula(model.formula),
                             data = data.frame(reftb[, !"i"], row.names = sample))



  ############################################
  # if data are batch corrected with sva
  if (!is.null(sva.obj)) {
    nds <- ncol(design.mat)
    nsv <- ncol(sva.obj$sv)
    svcol <- paste0("SV.", seq(nsv))
    design.mat <- cbind(design.mat, sva.obj$sv)
    colnames(design.mat)[(nds+1):(nds+nsv)] <- svcol
  }



  ############################################
  # clean column names in design matrix
  colnames(design.mat) <- gsub(paste0("^", fit.para[1]), "", colnames(design.mat))
  cols <- c(sample, fit.para, variable, value)
  dtl <- data[, ..cols] %>% unique %>% setorderv(., fit.para)
  dtw <- dcast.data.table(dtl, get(variable) ~ get(sample), value.var = value)
  mat <- as.matrix(dtw[, unlist(unique(dtl[, ..sample])), with = FALSE])
  rownames(mat) <- dtw[, variable]


  ############################################
  # lmFit procedure, w/out repeated measures
  #
  # Limma has a built-in approach for analyzing repeated measures data using duplicateCorrelation().
  # The model can handle a single random effect, and forces the magnitude of the random effect to
  # be the same across all genes.
  # Estimate linear mixed model with a single variance component
  # Fit the model for each gene,
  if (repeated.measures) {
    if (is.null(blocking.parameter)) stop("Blocking parameter not provided!")
    dupcor <- duplicateCorrelation(mat, design.mat,
                                   block = unlist(reftb[, ..blocking.parameter], use.names = FALSE))

    # convert refer3ence table to dataframe
    # rownames must match mat colnames
    # reftb %<>% .[, group := do.call(paste, c(.SD, sep = "___")), .SDcols = fit.para] %>%  data.frame(., row.names = "group")
    reftb %<>% data.frame(., row.names = sample)

    ## perform fitting
    fit0 <- lmFit(mat, design.mat,
                  block = unlist(reftb[, blocking.parameter]),
                  correlation = dupcor$consensus.correlation, na.action = na.exclude
    )
  } else {
    fit0 <- lmFit(mat, design.mat, na.action = "na.exclude")
  }


  ############################################
  # store fit residual
  fit_residuals <- residuals(fit0, mat) %>%
    as.data.table(., keep.rownames = variable) %>%
    melt(.,
         id.vars = variable,
         variable.name = sample,
         value.name = "Residual"
    ) %>%
    merge.data.table(., unique(dtl[, c(..sample, ..fit.para)]), by = sample)



  ############################################
  ## generate all possible pairwise combinations of condition levels
  if (pairwise.contrasts |  !is.null(manual.contrasts)) {
    para <- fit.para[1]
    para.lvl <- levels(unlist(dtl[, ..para]))
    if (uniqueN(para.lvl) == 1) {
      MSqb2:::.logg(ERROR, glue("Parameter {para} has 1 level in the design matrix. No contrast can be defined!"))
    } else {
      prs <- combn(para.lvl, 2)
    }


    # following condition levels better to be always in the denominator:
    # pairwise.denominator = c("wt", "wildtype", "wildt-ype", "cntrl", "control", "ctl", "untreated"),
    prs %<>% apply(., 2, \(x) {
      if (toupper(x[1]) %in% toupper(pairwise.denominator))
        replace(x, c(1, 2), x[c(2, 1)]) else x})


    # contrasts can be reversed if needed!
    if (reverse.contrasts) prs <- as.matrix(prs[2:1, ])
    comp <- apply(prs, 2, paste, collapse = " - ")

    # all possible contrasts
    all.contrasts <- makeContrasts(levels = design.mat, contrasts = comp)
    if (!is.null(sva.obj)) {
      all.contrasts <- rbind(all.contrasts,
                             matrix(
                               rep(
                                 rep(0, ncol(fit0$coefficients) - nrow(all.contrasts)),
                                 each = ncol(all.contrasts)
                               ),
                               ncol = ncol(all.contrasts))
      )
    }
  }



  ############################################
  # if all pairwise contrasts to be calculated
  if (pairwise.contrasts) cnts <- all.contrasts



  ############################################
  # choose pairwise contrasts
  if (choose.contrasts) {
    comp.p <- select.list(colnames(all.contrasts),
                          graphics = FALSE, multiple = TRUE,
                          title = "The following pairwise contrasts are available. Select the one(s) you are interested in."
    )
    cnts <- makeContrasts(levels = design.mat, contrasts = comp.p)
  }



  ############################################
  # generate manual pairwise contrasts
  if (!is.null(manual.contrasts)) { # handle the situation in which some manual.contrastsrasts have names and some don't
    manual.contrasts <- .rnManCont(manual.contrasts) # check names
    mcnt <- makeContrasts(levels = design.mat, contrasts = manual.contrasts)
    colnames(mcnt) <- names(manual.contrasts)
    if (!is.null(cnts)) {
      if (pairwise.contrasts) { # check if there is overlap between pairwise amd manual contrasts
        ix <- NULL
        for (i in seq(ncol(mcnt))) {
          ix <- c(
            ix,
            apply(cnts, 2, function(x) x == mcnt[, i]) %>%
              apply(., 2, all) %>%
              which()
          )
        }
        # remove overlapped contrasts from pairwise set
        if (length(ix) > 0) cnts <- cnts[, -ix]
      }

      cnts <- merge(cnts, mcnt, by = 0)
      rownames(cnts) <- cnts$Row.names
      cnts <- cnts[, setdiff(colnames(cnts), "Row.names")]
    } else {
      cnts <- mcnt
    }
  }



  ############################################
  # generate complex contrasts
  if (make.complex.contrasts) {
    message(
      "Given A, B, C and D to be four condition levels, \'B - C\' is a simple pairwise contrast.\n",
      "A complex contrast, on the other hand, can be defined like this: \n",
      "0.5*(A+B) - 0.5(C+D) i.e. contrast of average values\n",
      "or this\n",
      "0.5*(A-B) - 0.5(C-D) or the so-called delta-delta contrast.\n",
      "In the following you can choose the condition levels to be on the left side of the contrast ",
      "equation, like A and B in the example above, and those to be on the right side, like C and D in the example. ",
      "You must also determine if you are interested in a 'delta-delta' or an 'average' contrast!\n",
      "The complete equation will be printed on the console."
    )

    mkCont <- 1
    complex.contrasts <- NULL
    while (mkCont == 1) {
      cc <- select.list(c("delta-delta", "average"),
                         graphics = FALSE, multiple = FALSE,
                         title = "Do you want a 'delta-delta' (i.e. A-B vs C-D) \nor an 'average' (i.e. mean(A+B+...) vs mean(C+D+...)) contrast?"
      )

      if (cc == "delta-delta") {

        cat("Considerring the formula (A-B vs C-D), ...")
        ll1 <- select.list(colnames(design.mat),
                          graphics = FALSE, multiple = FALSE,
                          title = "choose Condition level to be passed to A.")
        ll2 <- select.list(setdiff(colnames(design.mat), ll1),
                          graphics = FALSE, multiple = FALSE,
                          title = "choose Condition level to be passed to B.")
        rr1 <- select.list(setdiff(colnames(design.mat), c(ll1, ll2)),
                           graphics = FALSE, multiple = FALSE,
                           title = "choose Condition level to be passed to C.")
        rr2 <- select.list(setdiff(colnames(design.mat), c(ll1, ll2, rr1)),
                           graphics = FALSE, multiple = FALSE,
                           title = "choose Condition level to be passed to D.")
        ll <- c(ll1, ll2)
        rr <- c(rr1, rr2)

      } else { # case of "average" contrast

        cat("Considerring the formula mean(A+B+...) vs mean(C+D+...),")
        ll <- select.list(colnames(design.mat),
                          graphics = FALSE, multiple = TRUE,
                          title = "Choose Condition level to be passed to the LEFT side of the contrast equation.")
        rr <- select.list(setdiff(colnames(design.mat), ll),
                          graphics = FALSE, multiple = TRUE,
                          title = "Choose Condition level to be passed to the right side of the contrast equation.")
      }



      # only with equal weights (at the moment)!
      if (cc == "average") {
        comp.c <- paste(
          length(ll)^-1 %>% paste0(., "*(", paste(ll, collapse = "+"), ")"),
          length(rr)^-1 %>% paste0(., "*(", paste(rr, collapse = "+"), ")"),
          sep = " - "
        )
      } else {
        comp.c <- paste(
          paste0("(", paste(ll, collapse = "-"), ")"),
          paste0("(", paste(rr, collapse = "-"), ")"),
          sep = " - "
        )
      }


      cat("The complex contrast formula is:\n", comp.c)
      complex.contrasts <- cbind(
        complex.contrasts,
        makeContrasts(levels = design.mat, contrasts = comp.c)
      )

      mkCont <- menu(c("yes", "no"), title = "Do you want to make another complex contrast?")
    }

    if (!is.null(cnts)) {
      cnts <- t(rbind(t(cnts), t(complex.contrasts)) )
    } else {
      cnts <- complex.contrasts
    }
  }



  if (add.vs.inContrast) colnames(cnts) <- gsub(" - ", " vs. ", colnames(cnts))


  ############################################
  # contrast fit and apply ebayes
  if (!time.course) fit_eb <- contrasts.fit(fit0, as.matrix(cnts)) else fit_eb <- fit0

  if (!is.logical(limma.trend)) {
    if (!limma.trend %in% names(data)) {
      MSqb2:::.logg(level = "WARN", msg = "limma.trend must either be logical or either of 'PSM.count' and 'PSM.mean'. By default it will be set to TRUE")
      limma.trend <- TRUE
    } else {
      PSM.df <- unique(data[, c(variable, limma.trend), with = FALSE]) %>% data.frame(., row.names = variable)
      limma.trend <- PSM.df[rownames(fit_eb$coefficients), limma.trend] %>% log2()
    }
  }
  fit_eb <- eBayes(fit_eb, trend = limma.trend, robust = limma.robust)


  if (removeTotalMissing) {
    ix <- which(apply(fit_eb$coefficients, 1, function(x) sum(is.na(x))) != ncol(fit_eb$coefficients))
    fit_eb <- fit_eb[ix, ]
    fit0 <- fit0[ix, ]
  }







  message("Current design matrix and contrasts ...")
  print(design.mat)
  if (!is.null(cnts)) print(cnts)

  ######################## from LIMMA guide
  ## Array Quality Weights
  # arrayw <- arrayWeights(mat)
  # barplot(arrayw, xlab="Array", ylab="Weight", col="white", las=2)
  # fitw <- lmFit(MAlms, design.mat, weights=arrayw)
  # fit2w <- contrasts.fit(fitw, contrast)
  # fit2w <- eBayes(fit2w)
  # topTable(fit2w, coef=1, number = 20, adjust="BH")


  # list2env(list("fit0" = fit0, "contrast" = contrast, "design.mat" = design.mat), envir = .GlobalEnv)
  save(fit_eb, fit0, fit_residuals, file = file.path(path, "fit_object.RData"))

  fl <- file.path(path, "fit_model.txt")
  sink(fl)
  cat("Design formula:\n==========\n")
  cat(model.formula)
  cat("\n\nDesign matrix:\n==========\n")
  print(design.mat)
  if (!is.null(cnts)) {
    cat("\nContrasts:\n==========\n")
    print(cnts)
  }
  sink()


  if (!is.null(cnts)) {
    write.table(cnts, row.names = FALSE, sep = "\t",
                file.path(path, paste0(prefix, "contrasts", suffix, ".txt"))
    )
  }

  options(warn = oldw)

  list(
    "fit_eb" = fit_eb,
    "fit_residuals" = fit_residuals,
    "fit0" = fit0
  ) %>% return(.)
}


## check manual contrasts names
.rnManCont <- function(mcx) {
  if (!is.list(mcx)) mcx <- as.list(mcx)
  if (is.null(names(mcx))) names(mcx) <- mcx
  nx <- which(names(mcx) == "")
  names(mcx)[nx] <- as.character(mcx[nx])
  return(mcx)
}
