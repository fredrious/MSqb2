
# %%%%%%%%%%%%%%%%%%%%%%%%
# import all objects from parent environment (parent env.)
.impTopEnv <- function() {
   list2env(mget(names(parent.frame(n = 2)), envir = parent.frame(n = 2)),
            envir = parent.frame(n = 1))
}


# %%%%%%%%%%%%%%%%%%%%%%%%
# export all objects to parent environment (parent env.)
.expTopEnv <- function() {
   list2env(mget(names(parent.frame(n = 1)), envir = parent.frame(n = 1)),
            envir = parent.frame(n = 2))
}

# %%%%%%%%%%%%%%%%%%%%%%%%
# import (expand) all variables in a list to the current environment
.expandList <- function(x) {
   list2env(mget(names(x), envir = as.environment(x)),
            envir = parent.frame())
}


# %%%%%%%%%%%%%%%%%%%%%%%%
# Assign variable in Parent Environment of a function
# either one step backward (n=2 in `%<-1%`) or 2 steps (n=3 in  `%<-2%`)
`%<-1%` <- function(x, y) { assign(x, y, pos = parent.frame(n = 2)) }
`%<-g%` <- function(x, y) { assign(x, y, pos = .GlobalEnv) }



# %%%%%%%%%%%%%%%%%%%%%%%%
## check build arguments
.checkBuildArgs <- function() {

   .impTopEnv() #import all objects from 1 env backward (parent.frame)

   MSqb2:::.loggAssert(assertCharacter(measurements.file))
   MSqb2:::.loggAssert(assertCharacter(metadata.file, len = 1))

   assert(checkCharacter(measurements.file.sheet, len = 1),
          checkInt(measurements.file.sheet, null.ok = TRUE),
          combine = "or")

   assert(checkCharacter(metadata.file.sheet, len = 1),
          checkInt(metadata.file.sheet, null.ok = TRUE),
          combine = "or")


   assert(checkChoice(ms.software, c("MQ", "PD")))

   MSqb2:::.loggAssert(assertCharacter(analysis.name, len = 1, null.ok = TRUE))
   MSqb2:::.loggAssert(assertCharacter(ms.software, len = 1, null.ok = FALSE))
   MSqb2:::.loggAssert(assertCharacter(prefix, len = 1, null.ok = TRUE))
   MSqb2:::.loggAssert(assertCharacter(suffix, len = 1, null.ok = TRUE))

   MSqb2:::.loggAssert(assertLogical(add.date.tag))
}




# %%%%%%%%%%%%%%%%%%%%%%%%
## logging wrapper file (issue with colored text in appender_file. )
.logg <- function(level, msg) {

   # define/call log file
   log_file <-
   if ("log_file" %in% ls(envir = .GlobalEnv)) {
      get("log_file", envir = .GlobalEnv)
   } else {
      "_temp_log_file.log"
   }

   # get level
   l.level <-  substitute(level) %>% as.character %>% tolower()
   level <- l.level %>% toupper

   # include all logging levels
   log_threshold(TRACE)

   #!! logging with color-coding does not work properly with append_file/tee!
   # I defined therefore two separate logging layout for console and file.
   # print in file:
   logger.file <-
      layout_glue_generator(
         format = '---\n••• MSqb {level}:\n{msg}\n')
   log_appender(appender_file(log_file,))
   log_layout(logger.file)
   do.call(glue("log_{l.level}"), list(msg))

   # print in console:
   logger.console <-
      layout_glue_generator(
         format = '{colorize_by_log_level(paste("---\n••• MSqb", level, ":\n", msg), levelr)}\n')
   log_appender(appender_console)
   log_layout(logger.console)
   do.call(glue("log_{l.level}"), list(msg))
   # options(show.error.messages = FALSE)
   if (tolower(level) %in% c("fatal", "error")) stop(call. = FALSE)
}




# %%%%%%%%%%%%%%%%%%%%%%%%
## make sub directories
.mk.dir <- function(path, sub.dir, prefix = NULL, analysis.name = NULL,
                    suffix = NULL, add.date.tag = NULL) {
   for (i in seq_along(sub.dir)) {
      sdir.p <- file.path(path,
                          paste(c(prefix, sub.dir[i], analysis.name, suffix, add.date.tag),
                                collapse = "_"))
      ifelse(!dir.exists(sdir.p), dir.create(sdir.p, recursive = TRUE), FALSE)
      paste0(sub.dir[i], ".path") %<-1% sdir.p
   }
}




# %%%%%%%%%%%%%%%%%%%%%%%%
## find input file (used inside .check.dir)
.find.file <- function(input, output, toData, where) {

   # input %<>% .[setdiff(names(input), names(output))] ## check this !!!

   fls <- lapply(input, function(x) grep(pattern = x,
                                  list.files(where,
                                             all.files = TRUE,
                                             recursive = TRUE,
                                             full.names = TRUE),
                                  value = TRUE))
   fls <- fls[which(!sapply(fls, function(x) identical(x, character(0))))]
   mx <- which(sapply(fls, length) > 1)
   if (length(mx) > 0) {
      MSqb2:::.logg(level = FATAL,
            glue("Multiple matches were found for the following input file(s):\n",
                 "* {glue_collapse(unlist(fls[mx]), sep = '\n* ')}",
                 "\n\nPlease provide the complete file path."))
   }

   if (length(fls) > 0) {
   list2env(list("output" = append(output, fls),
                 "input" = input,
                 "toData" = toData),
            envir = parent.frame())
   }
}



# %%%%%%%%%%%%%%%%%%%%%%%%
## check input data directories and copy to Data if approved by user
.check.file <- function(input, dpth = Data.path, toData = FALSE, interactive = TRUE) {

   input %<>% lapply(., function(x) !is.null(x)) %>% unlist %>% input[.]
   output <- list()
   ninp <- length(input)


   ## check if input files are full paths and exist
   idx <- which(file.exists(unlist(input)))
   if (length(idx) > 0) {
      output[names(input)[idx]] <- input[idx]
      input <- input[setdiff(names(input), names(output))]
      toData <- TRUE
   }


   ## check if input files exist in dpth path
   if (length(input) > 0) .find.file(input, output, toData, dpth)


   ## check if input files not exist in dpth path, but do exist in working directory
   if (length(output) != ninp) .find.file(input, output, toData, getwd())


   ## stop if files could not be found!
   if  (length(output) != ninp) {
      nafile <- setdiff(names(input), names(output)) # the file that was not found!

      if (basename(dpth)  != "Scripts") {
      MSqb2:::.logg(level = FATAL,
            glue("Following file(s) could not be found in the working directory.\n",
                 "{glue_collapse(nafile, sep = '\n')}\n",
                 "Either copy them in the working directory (preferably in 'Data' directory) ",
                 "or provide the full path.") )
      } else {
         MSqb2:::.logg(level = WARN,
               glue("The following config file does not exist or could not be found in the working ",
                    "directory.\n {glue_collapse(input, sep = '\n')}\n",
                    "By default, a config file will be generated using the internal config file template. ",
                    "Alternatively you can provide the correct path and rerun msqb_config.") )
         output <- ""
      }
      Sys.sleep(0.1)
   }


   if (toData) {
      idx <- which(lapply(output, dirname) != dpth)
      if (length(idx) > 0 & interactive) {
         MSqb2:::.logg(TRACE, glue(
            "Input files are recommended (but not required!) to be stored ",
            "in the {basename(dpth)} directory: \n{dpth} \n\n",
            "Do you want to copy the following files in the 'Data' directory?\n\n",
            paste(output[names(idx)], collapse = "\n")))

         cp.data <- select.list(c("yes", "no"))

         if (cp.data == "yes") {
            lapply(output[names(idx)], function(x) file.copy(from = x, to = dpth, overwrite = FALSE))
            MSqb2:::.logg(SUCCESS, glue(
               "Following file(s) were copied to '{basename(dpth)}' sub-directory:\n",
               "{glue_collapse(output[names(idx)], sep = '\n')}\n",
               "\nThese files will be utilised throughout the workflow and in the future analysis.\n"))

            output[names(idx)] <- lapply(sapply(output[names(idx)], basename), function(x) file.path(dpth, x))
         }
      }
   }
   return(output)
}




# %%%%%%%%%%%%%%%%%%%%%%%%
read.file <- function(file, sheet = 1) {

   ext <- toupper(strsplit2(file, "\\.")[-1])
   if (any(ext %in% c("XLSX", "XLS"))) {
      dt <- readxl::read_excel(file, sheet = sheet) %>% as.data.table(.)
   } else {
      dt <- fread(file, header = TRUE, sep = "\t", stringsAsFactors = TRUE)
   }
   MSqb2::char2fact(dt) %>% 
   return()
}




# %%%%%%%%%%%%%%%%%%%%%%%%
## import config parameters from file (call from msqb_config)
.importConfigfile <- function(conf.fl, Scripts.path, analysis.name,
                              add.date.tag, whichConf, conf.args) {

   if (!is.null(conf.fl)) conf.pth <- .check.file(conf.fl, dpth = Scripts.path, toData = FALSE)
   # if (exists("conf.pth") && conf.pth != "") conf.pth <- basename(conf.fl)
   if ((!exists("conf.fl") || is.null(conf.fl) ||
        (exists("conf.pth") && conf.pth == "" ) )) { #if configuration file does not exist, create one!
      conf.pth <- file.path(
         Scripts.path,
         paste0(
            ifelse(whichConf == "workflow", "msqb_config", "msqb_config_viz"),
            ifelse(is.null(analysis.name), "", paste0("_", analysis.name)),
            add.date.tag, ".R") )
      Sys.sleep(0.2)
      makeConfigFile(conf.file = conf.pth, whichConf = whichConf)
   }

   suppressWarnings(
      if (length(within(conf.args, rm(interactive))) > 0) {
         MSqb2:::.logg(TRACE,
              glue(
            "Parameters manually passed to the msqb_config function will substitude ",
            "the parameters in the config.file (if provided). These are as follows:\n",
            "{glue_collapse(names(within(conf.args, rm(interactive))), sep = '\n')} \n"))
      }
   )

   MSqb2:::.logg(INFO, glue("Config parameters are stored under: {conf.pth}"))
   conf.pth <- unlist(conf.pth)
   eval(parse(text = readLines(conf.pth)), envir = parent.frame())
   return(conf.pth)
}




# %%%%%%%%%%%%%%%%%%%%%%%%
## combine logger and checkmate, so the assertion message will have the default msqb format
.loggAssert <- function(xpr) {
   .impTopEnv()
   .collAssert = makeAssertCollection()
   xpr <- deparse(substitute(xpr)) %>% gsub(" ", "", .) %>% paste(., collapse = "")
   xpr <- paste0(substr(xpr, 1, nchar(xpr)-1), ",add=.collAssert)")
   eval(parse(text = xpr))
   if (!.collAssert$isEmpty()) {
      MSqb2:::.logg(error, skip_formatter(.collAssert$getMessages()))
   }
}




# %%%%%%%%%%%%%%%%%%%%%%%%
## convert factor columns to character columns in data.table
fact2char <- function(dt) {
  changeCols <- c(names(Filter(is.factor, dt)))
  if (length(changeCols) > 0) {
    dt[, (changeCols) := lapply(.SD, as.character), .SDcols = changeCols]
  }
  return(dt)
}




# %%%%%%%%%%%%%%%%%%%%%%%%
#' character to factor converter
#'
#' Converts data.frame or data.table columns of class character to factor.
#'
#' @import data.table
#' @param dt Input data of class data.table or data.frame.
#' @return The input data \code{dt} with columns of class character converted to factor.
#' @seealso \code{\link{int2fact}}
#' @examples
#' DT <- data.table(x = c(1., 2., 3.), y = c(4L, 5L, 6L), z = c("a", "b", "NA"))
#' # DT <- data.frame(x = c(1.,2.), y = c(1L,2L), z = c("1","2"))
#' sapply(DT, class)
#' MSqb2::char2fact(DT)
#' sapply(DT, class)
char2fact <- function(dt) {
  if (!any(class(dt) %in% c("data.table", "data.frame"))) {
    stop("Input must be of class data.table or data.frame.")
  } else if (!is.data.table(dt)) {
    rn <- rownames(dt)
    setDT(dt)
    setback2df <- TRUE
  } else {
    setback2df <- FALSE
  }

  changeCols <- names(Filter(is.character, dt))
  if (length(changeCols) > 0) {
    dt[, (changeCols) := lapply(.SD, \(x) {
      is.na(x) <- x %in% c("NA", "<NA>")
      x <- factor(x, levels = unique(x)) 
    }), .SDcols = changeCols]
  }

  # drop excluded levels
  fc <- names(Filter(is.factor, dt))
  dt[, (fc) := lapply(.SD, droplevels), .SDcols = fc]
  # sort based on alphanumeric -> natural (lexicographic)
  dt[, (fc) := lapply(.SD, \(x) factor(x, levels = unique(x) %>% stringr::str_sort(., numeric = TRUE) )), .SDcols = fc]

  
  if (setback2df) dt <- data.frame(dt, row.names = rn)
  return(dt)
}




# %%%%%%%%%%%%%%%%%%%%%%%%
#' integer to factor converter
#'
#' Converts data.frame or data.table columns of class integer to factor.
#'
#' @import data.table
#' @param dt Input data of class data.table or data.frame.
#' @return The input data \code{dt} with columns of class integer converted to factor.
#' @seealso \code{\link{MSqb2::char2fact}}
#' @examples
#' DT <- data.table(x = c(1., 2.), y = c(1L, 2L), z = c("1", "2"))
#' # DT <- data.frame(x = c(1.,2.), y = c(1L,2L), z = c("1","2"))
#' sapply(DT, class)
#' int2fact(DT)
#' sapply(DT, class)
int2fact <- function(dt) {
  if (!any(class(dt) %in% c("data.table", "data.frame"))) {
    stop("Input must be of class data.table or data.frame.")
  } else if (!is.data.table(dt)) {
    rn <- rownames(dt)
    setDT(dt)
    setback2df <- TRUE
  } else {
    setback2df <- FALSE
  }

  changeCols <- c(names(Filter(is.integer, dt)))
  if (length(changeCols) > 0) {
    dt[, (changeCols) := lapply(.SD, as.factor), .SDcols = changeCols]
  }

  if (setback2df) dt <- data.frame(dt, row.names = rn)
  return(dt)
}





# %%%%%%%%%%%%%%%%%%%%%%%%
## match model formula with experimental design
.match_ModelFormula_metadata <- function(dsgn, frm) {
   fit.para <- c(all.vars(as.formula(frm)))
   if (any(!fit.para %in% names(dsgn))) {
      xpara <- setdiff(fit.para, names(dsgn))
      MSqb2:::.logg(TRACE,
         glue(
            "The following parameter(s) that used in the model formula do ",
            "not match the column names in the metadata! \n\n",
            "missing parameter(s): {xpara}",
            "\nmodel formula: {frm}",
            "\n\nDo you want to stop the process and modify the metadata ",
            "or model formula? By selecting 'No' the missing parameter(s) ",
            "will be removed from  model formula and the workflow continues. ",
            "Note that for complex formulae this might result in a false formula. ",
            "It is therefore recommended to double-check the model formula or rerun ",
            "the analysis after correcting the input data or the model formula."
         ))

      mtch <- select.list(c("yes", "no"))

      if (mtch == "yes") {
         MSqb2:::.logg(ERROR, "Error in matching model formula and cathegorical variables in the metadata.")
      } else frm <- gsub(xpara, "NULL", frm)
   }
   return(frm)
}





# %%%%%%%%%%%%%%%%%%%%%%%%
## import config parameters from file
readConfigPara <- function(build.para, config.para, config.file, config.type, ...) {
   if (!exists("build.para", mode = "list")) {
      MSqb2:::.logg(FATAL, glue(
         "build.para not provided. ",
         "Please run msqb_build() to generate build parameters (see vignette)."
      ))
   } else .expandList(build.para)

   conf.str <- switch(config.type,
                      "viz" = c("config.vis.file", "config.para.viz"),
                      "wf" = c("config.file", "config.para"))


   if ( (exists("config.file") && !is.null(config.file)) &
        (exists("config.para") && !is.null(config.para)) ) {
      MSqb2:::.logg(TRACE, glue(
         "Both {conf.str[1]} and {conf.str[2]} have been provided! ",
         "Do you want the parameters stored in {conf.str[1]} to be loaded? ",
         "By chosing NO the parameters in {conf.str[2]} will be loaded and ",
         "a new {conf.str[1]} will be generated."
      ))
      conf.res <- select.list(c("yes", "no"))
      if (conf.res == "yes") config.para <- NULL else config.file <- NULL
   }


   ## either read or generate config parameters
   conf.args <- as.list(substitute(list(...)))[-1L]
   if (!exists("config.para", mode = "list") || is.null(config.para) || length(conf.args) > 0) {
      conf.para <-
         switch(
            config.type,
              "wf"  = msqb_config(build.para = build.para, config.file = config.file, ...),
              "viz" = msqb_config_viz(build.para = build.para, config.viz.file = config.file, ...)
         )
   } else conf.para <- config.para

   # expand conf.para
   .expandList(conf.para)
   # remove junk
   suppressWarnings(rm(conf.para, conf.str, config.file, config.type, conf.res, config.para))
   # export all to top env.
   .expTopEnv()
}




# %%%%%%%%%%%%%%%%%%%%%%%%
## Add BioRep if not already in the metadata
.check.pheno <- function(pdt) {

   ## --- paired design:
   ## Paired samples occur when we compare two treatments and each sample given one treatment
   ## is naturally paired with a particular sample given the other treatment.
   ## (paired samples with blocks of size 2, 3, ...)

   ## --- factorial design
   ## Factorial designs are those where more than one experimental dimension is being varied
   ## and each combination of treatment conditions is observed.
   ## -> here make a new condition column by combining (pasting) the dimensions together.

   ## limma 9.7 Multi-level Experiments
   ## the case with duplicateCorrelation


   if (!"BioRep" %in% names(pdt)) { # if no BioRep in data
      if (!"SampleID" %in% names(pdt)) { # if no BioRep AND sampleID in data
         MSqb2:::.logg(WARN, glue("There are no 'SampleID' and 'BioRep' columns in the metadata. ",
                          "By default, each raw will be considered as a biological replicate."))
         pdt[, SampleID := sprintf("sample.%0*s", floor(log10(max(.I)))+1, .I) ]
      } else {
         MSqb2:::.logg(WARN, glue("There is no 'BioRep' column in the metadata. ",
                          "By default, each sample will be considered as a biological replicate."))
      }
      pdt[, BioRep := sprintf("BioRep.%0*s", floor(log10(max(.I)))+1, .GRP), by = SampleID]
   }








   if (any(pdt[, .(ConditionLvls = uniqueN(Condition)), by = SampleID][, 2] > 1)) {

      MSqb2:::.logg(INFO, glue(
         "Multiple condition levels found for the same samples. A 'repeated measures' ",
         "design will be considered for the statistical analysis." ))

      SmplConds <- pdt[, .(ConditionLvls = toString(Condition)), by = SampleID]
      print(SmplConds)
      # also add to log file
      sink(log_file, append = TRUE)
      print(SmplConds)
      sink()

      repeated.measures <- TRUE
      blocking.para <- "SampleID"
   }




   if (!"BioRep" %in% names(pdt) & !"TechRep" %in% names(pdt)) {
      MSqb2:::.logg(WARN, glue("There are no 'BioRep' and 'TechRep' columns in the metadata. ",
                       "By default, each sample will be considered as a biological replicate."))
      pdt[, BioRep := paste0("BioRep", seq_len(.N)), by = Condition]
   }

   if (!"BioRep" %in% names(pdt) & "TechRep" %in% names(pdt)) {
      MSqb2:::.logg(WARN, glue("There is no 'BioRep' column in the metadata. ",
                       "By default, 'BioRep' column will be created by sequencing along ",
                       "'Condition' levels grouped by 'TechRep' levels."))
      pdt[, BioRep := paste0("BioRep", seq_len(.N)), by = Condition]
      pdt[, BioRep := paste0("TechRep", seq_along(Condition)), by = "TechRep"]
   }

   MSqb2:::.logg(INFO, glue("metadata:\n"))
   print(pdt)
   # also add to log file
   sink(log_file, append = TRUE)
   print(pdt)
   sink()

   return(invisible(pdt))
}




# %%%%%%%%%%%%%%%%%%%%%%%%
## wide to long format & matrix to long
.w2l <- function(dw,
                 val = "Abundance",
                 col.p,
                 row.p,
                 sp = "&.&.&") {
  if (is.matrix(dw)) {
    dw <- as.data.table(dw, keep.rownames = "nm")
    dw[, (row.p) := tstrsplit(nm, sp, fixed = TRUE)]
    dw$nm <- NULL
  }

  dl <- melt(dw,
             id.vars = row.p,
             variable.name = "dummy",
             value.name = val
  )

  dl[, (col.p) := tstrsplit(dummy, sp, fixed = TRUE)]
  dl[, dummy := NULL]

  dl <- MSqb2::char2fact(dl)
  return(dl)
}



# %%%%%%%%%%%%%%%%%%%%%%%%
## long to wide format & long to matrix
.l2w <- function(dl,
                 col.p,
                 row.p,
                 val = "Abundance",
                 sp = "&.&.&",
                 asMat = FALSE,
                 Mat.rownm = NULL, ...) { # ... for fun.aggregate

  # define casting formula
  form <- paste0(
    ifelse(length(row.p) > 1, paste(row.p, collapse = "+"), row.p),
    "~",
    ifelse(length(col.p) > 1, paste(col.p, collapse = "+"), col.p)
  )

  dw <- dcast.data.table(dl,
                         formula = form,
                         value.var = eval(val),
                         sep = sp, ...
  ) # cast input data -- wide format

  if (!is.null(Mat.rownm)) {
    rownm <- dw[, do.call(paste, c(.SD, sep = sp)), .SDcols = Mat.rownm]
  }
  if (asMat) {
    dm <- as.matrix(dw[, !..row.p], rownames = rownm)
  } else {
    dm <- dw
  }
  return(dm)
}






# %%%%%%%%%%%%%%%%%%%%%%%%
## !!! FOLLOWING SYMBOLS ARE NOT ALLOWED IN COLUMN NAMES AND INSIDE THE COLUMNS.
##     This is to avoid any confusion with the symbols used in the design formula.

cleanDataPara <- function(ddt, complete.check = FALSE, # export2env = FALSE
                           subsymDT = data.table(
                             symout = c("\\~", "\\+", "\\:", "\\|", "\\-", "\\(", "\\)", "\\/", " "),
                             symin = c("", ".", "", "", "_", "", "", ".", "")
                           )) {

  if (prod((c("symout", "symin") %in% names(subsymDT))) != 1) {
    subsymDT[, var := c("symout", "symin")] %>% meltsub.f(.) -> subsymDT
  }



  .meltsub.f <- function(dt) {
    return(
      melt(dt, id.vars = "var")[, -"variable"] %>%
        .[, ii := seq_along(value), by = "var"] %>%
        dcast(., formula = ii ~ var) %>%
        .[, -"ii"]
    )
  }



  if (is.data.table(ddt)) {
    nmDT <- data.table(original.names = names(ddt)) # old names with punctuations
    for (is in seq(nrow(subsymDT))) {
      names(ddt) <- gsub(subsymDT$symout[is], subsymDT$symin[is], names(ddt))
    }
    names(ddt) <- make.names(names(ddt))
    names(ddt) <- sapply(names(ddt), function(x) {
      strsplit2(x, "\\.") %>%
        .[nchar(.) > 0] %>%
        paste(., collapse = ".")
    })
    nmDT[, new.names := names(ddt)] # new names with punctuations

    if (complete.check) {
      for (is in seq(nrow(subsymDT))) {
        for (para in nmDT$new.names) {
          ddt[, (para) := lapply(.SD, function(x) gsub(subsymDT$symout[is], subsymDT$symin[is], x, perl = TRUE)), .SDcols = para]
        }
      }
    }

    # message("\n\nSpecial charachters in the column names of 'Annotation.dt' were substituted according to 'subsymDT' table.")
    # print(subsymDT[, .(symout, symin)])
    # message("!!! Parameters used in design formula for test statistics must match the new names.")
    # message("Check the original and new names in the table below.")
    # print(nmDT)

    # if (export2env) {
    #   list2env(list("subsymDT" = subsymDT,
    #                 "ref.name.tbl" = nmDT),
    #          envir = .GlobalEnv)
    # }
  }

  if (is.vector(ddt)) {
    for (is in seq(nrow(subsymDT))) {
      ddt <- gsub(subsymDT$symout[is], subsymDT$symin[is], ddt)
      nmDT <- NULL
    }
  }

  dout <- list("cln.nm.dt" = ddt, "ref.name.tbl" = nmDT)

  return(dout)
}





# %%%%%%%%%%%%%%%%%%%%%%%%
## check color palettes, install if needed

# test: plts <- list(Condition = 'Set1', BioRep = 'Darjeeling2', sd = "Heat", hh = "npg")
.check.palettes <- function(plts, return.pkg = FALSE) {

  if (is.list(plts)) plts %<>% unlist(.)

  plt.l <- list(
    RColorBrewer = c('BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn',
                     'Spectral', 'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1',
                     'Set2', 'Set3', 'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys',
                     'Oranges', 'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds',
                     'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd'),

    viridis = c('cividis', 'inferno', 'magma', 'mako', 'plasma', 
                'rocket', 'turbo', 'unemp', 'viridis'),

    wesanderson = c('BottleRocket1', 'BottleRocket2', 'Rushmore1', 'Rushmore', 'Royal1',
                    'Royal2', 'Zissou1', 'Darjeeling1', 'Darjeeling2', 'Chevalier1',
                    'FantasticFox1', 'Moonrise1', 'IsleofDogs2', 'Moonrise2', 'Moonrise3',
                    'Cavalcanti1', 'GrandBudapest1', 'GrandBudapest2', 'IsleofDogs1'),

    ggpubr = c('npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', 'ucscgb', 'd3', 'locuszoom',
               'igv', 'uchicago', 'startrek', 'tron', 'futurama', 'rickandmorty', 'simpsons'),

    grDevices = c('Pastel 1', 'Dark 2', 'Dark 3', 'Set 2', 'Set 3', 'Warm', 'Cold', 'Harmonic',
                  'Dynamic', 'Grays', 'Light Grays', 'Blues 2', 'Blues 3', 'Purples 2',
                  'Purples 3', 'Reds 2', 'Reds 3', 'Greens 2', 'Greens 3', 'Oslo', 'Purple-Blue',
                  'Red-Purple', 'Red-Blue', 'Purple-Orange', 'Purple-Yellow', 'Blue-Yellow',
                  'Green-Yellow', 'Red-Yellow', 'Heat', 'Heat 2', 'Terrain', 'Terrain 2',
                  'Viridis', 'Plasma', 'Inferno', 'Rocket', 'Mako', 'Dark Mint', 'Mint',
                  'BluGrn', 'Teal', 'TealGrn', 'Emrld', 'BluYl', 'ag_GrnYl', 'Peach', 'PinkYl',
                  'Burg', 'BurgYl', 'RedOr', 'OrYel', 'Purp', 'PurpOr', 'Sunset', 'Magenta',
                  'SunsetDark', 'ag_Sunset', 'BrwnYl', 'YlOrRd', 'YlOrBr', 'OrRd', 'Oranges',
                  'YlGn', 'YlGnBu', 'Reds', 'RdPu', 'PuRd', 'Purples', 'PuBuGn', 'PuBu', 'Greens',
                  'BuGn', 'GnBu', 'BuPu', 'Blues', 'Lajolla', 'Turku', 'Hawaii', 'Batlow',
                  'Blue-Red', 'Blue-Red 2', 'Blue-Red 3', 'Red-Green', 'Purple-Green', 'Purple-Brown',
                  'Green-Brown', 'Blue-Yellow 2', 'Blue-Yellow 3', 'Green-Orange', 'Cyan-Magenta',
                  'Tropic', 'Broc', 'Cork', 'Vik', 'Berlin', 'Lisbon', 'Tofino', 'ArmyRose',
                  'Earth', 'Fall', 'Geyser', 'TealRose', 'Temps', 'PuOr', 'RdBu', 'RdGy', 'PiYG',
                  'PRGn', 'BrBG', 'RdYlBu', 'RdYlGn', 'Spectral', 'Zissou 1', 'Cividis', 'Rom')
  )

  reserved.plts <- c('Set1', 'Spectral', 'Accent', 'Dark2',
                     'Paired','YlGnBu', 'RdYlBu', 'Set2', 'Set3')


  # check which palettes are not allowed
  if (any(!plts %in% unlist(plt.l))) {
    MSqb2:::.logg(ERROR,
                 glue("Color palette not allowed:\n",
                      "{glue_collapse(setdiff(plts, plt.l), sep = '\n')}\n"))
  }


  # check which packeges needed to be loaded/installed
  pltpkgs <- sapply(plts, \(x) grep(x, plt.l)[1] ) %>% plt.l[.] %>% names()
  ipk <- which(!pltpkgs %in% rownames(installed.packages()))

  # if ipk not empty, some packages are needed to be installed
  if (length(ipk) != 0) {

    instpkg <- pltpkgs[ipk]

    if (interactive) {
      MSqb2:::.logg(TRACE,
                   glue("Following color palettes for visualization are from package ",
                        "which are not installed:\n\n",
                        "{glue_collapse(plts[ipk], sep = '\n')}\n\n",
                        "required packages (respectively):\n\n",
                        "{glue_collapse(instpkg, sep = '\n')}\n\n",
                        "Do you want to install these packages?\n",
                        "By selecting 'no', these palettes will be replaced by random color palettes. ",
                        "Alternatively you can terminate the process and change the palettes."))

      ans <- select.list(c("yes", "no"))

      if (ans == "yes") {

        install.packages(instpkg)
        if (!all(instpkg %in% rownames(installed.packages()))) { # check again
          library(BiocManager, quietly = TRUE)
          BiocManager::install(instpkg)
        }
        sapply(instpkg, require, character.only = TRUE)
        MSqb2:::.logg(SUCCESS, glue("Package(s) were successfully installed and loaded."))

      } else {
        plts[ipk] <- sample(setdiff(reserved.plts, plts[-ipk]), 2)
      }
    }
  }

  if (return.pkg) {
    data.table(pkg = sapply(plts, \(x) grep(x, plt.l)[1] ) %>% plt.l[.] %>% names(),
               plt = plts) %>%
      return(.)
  } else return(as.list(plts))

}





# %%%%%%%%%%%%%%%%%%%%%%%%
## pick colors from palettes

colorpicker <- function(plt, n, ...) {
  oldw <- getOption("warn")
  options(warn = -1)

  # pltdt <- MSqb2:::.check.palettes(plt, return.pkg = TRUE)
  pltdt <- .check.palettes(plt, return.pkg = TRUE)
  plt <- pltdt[1, plt]
  pkg <- pltdt[1, pkg]

  if (pkg == "wesanderson") {

    cols <- as.vector(wesanderson::wes_palette(n = n, name = plt, type = c("continuous")))

  } else if (pkg == "ggpubr") {

    cols <- ggpubr::get_palette(k = n, palette = plt)

  } else if (pkg == "grDevices") {

    cols <- grDevices::hcl.colors(n = n, palette = plt, ...)

  } else if (pkg == "RColorBrewer") {

    getPalette = colorRampPalette(RColorBrewer::brewer.pal(n, plt))
    cols <- getPalette(n)
    # cols <- RColorBrewer::brewer.pal(n = n, name = plt)
    if (n < 3) cols <- cols[1:n]

  } else if (pkg == "viridis") {

    cols <- viridis::viridis_pal(option = tolower(plt))(n)

  }
  options(warn = oldw)

  return(cols)
}





# %%%%%%%%%%%%%%%%%%%%%%%%
## display global colors
.displayGlobalColors <- function(col.list, silent = FALSE, save.plot = FALSE, path = getwd()) {
  cols <- suppressWarnings(as.data.table(col.list))
  nms <- suppressWarnings(as.data.table(lapply(col.list, names)))
  .col2l <- function(d, ...) {
    return(d[, id := seq(nrow(d))] %>%
             melt.data.table(., id.vars = "id", ...) %>% .[, -"id"] %>% unique(.))
  }

  global.colors <-
    cbind(
      .col2l(cols, value.name = "cl"),
      .col2l(nms, value.name = "pr")[, -"variable"]
    ) %>%
    .[, y := seq_along(cl), by = variable] %>%
    ggplot(., aes(x = variable, y = y, col = cl, fill = cl, label = pr)) +
    geom_tile() +
    facet_wrap(. ~ variable, scales = "free") +
    geom_text(col = "black") +
    scale_color_identity() +
    scale_fill_identity() +
    theme_void()

  if (save.plot) {
    png(filename = file.path(path, "global_colors_scheme.png"))
    print(global.colors)
    dev.off()
  }
  if (!silent) {
    return(global.colors)
  }
}






# %%%%%%%%%%%%%%%%%%%%%%%%
## define facet formula for qc plots
.FacetFormula <- function(fct) {
  if (!is.list(fct)) stop("Input must be of class list, with 2 elements: 'rows' and 'cols'!")
  if (length(fct) > 2) stop("Input must be of class list, with 2 elements: 'rows' and 'cols'!")
  if (any(is.na(match(names(fct), c("rows", "cols"))))) {
    stop("Invalid element names! The name of the elements in the input list must be 'rows' and 'cols'")
  }

  elm <- setdiff(c("rows", "cols"), names(fct))
  if (length(elm) > 0) fct[[elm]] <- "."

  return(formula(paste(
    paste(fct$rows, collapse = "+"),
    paste(fct$cols, collapse = "+"),
    sep = "~"
  )))
}





# %%%%%%%%%%%%%%%%%%%%%%%%
## if feature.annotation.source is psm.file, call annotations from description column
.psmDescGenes <- function(x) {
  xx <- strsplit(x, "GN=")
  unlist(xx) %>% .[-1] %>% 
    sapply(., \(y) tstrsplit(y, " PE=", keep = 1)) %>% 
    unlist() %>% 
    paste(., collapse = ";")
}





# %%%%%%%%%%%%%%%%%%%%%%%%
## order data.table with alphanumeric columns
## idea from: https://stackoverflow.com/questions/49084625/sort-strings-by-numbers-inside-them
order.num <- function(dt, x) {
  xo <- dt[, ..x] %>% unlist()
  fct <- ifelse(is.factor(xo), TRUE, FALSE)
  if (fct) xo %>% unique() %>% as.character()
  ix <- xo %>% gsub("[^[:digit:]]", "", .) %>% as.numeric() %>% order() %>% xo[.]
  dt <- dt[order(match(get(x), ix))] 
  if (fct) dt[, (x) := factor(get(x), levels = ix)]
  return(dt)
}





# %%%%%%%%%%%%%%%%%%%%%%%%
## merge multi column data with lookup table
lookupSplit <- function(dpg, ref, input, output) {
  setkeyv(ref, cols = input)

  dpg[, rn := seq(nrow(dpg))] %>%
    # wide to long
    melt(., id.vars = "rn", variable.name = "inputGroup", value.name = input) %>%
    .[!is.na(get(input))] %>%
    # Set key for the join
    setkeyv(., input) %>%
    ref[.] %>%
    # wide to long
    dcast(., rn ~ inputGroup, value.var = output) %>%
    .[, rn := NULL] %>%
    return(.)
}





# %%%%%%%%%%%%%%%%%%%%%%%%
#' Subsets Top proteins/genes Based on Variance
#'
#' This function subsets the input data table based on the highest row-wise variance.
#' It can return either a specific number of top rows (`topN`) or a percentage
#' (`topNperc`) of the rows with the highest variance. Additional options allow
#' control over the minimum number of rows returned and the format of the output.
#' The combination of `row.p`, `col.p` and `val` will be used to reformat data via
#' function `dcast` from `data.table` package.
#'
#' @param dt Data table to be processed.
#' @param row.p Column name in `dt` representing the rows.
#' @param col.p Column name in `dt` representing the columns.
#' @param val Column name in `dt` representing the values.
#' @param topNperc Percentage of top variant rows to return. When specified, `topN` is ignored.
#' @param min500 Logical; if TRUE, a minimum of 500 rows are returned.
#' @param topN Integer; the number of top variant rows to return. When specified, `topNperc` is ignored.
#' @param complete.rows Logical; if TRUE, rows with missing values are omitted.
#' @param sep.col String; separator used in the matrix conversion.
#' @param returnLongFormat Logical; if TRUE, returns the data in long format.
#' @return A subset of `dt` with the top variable rows based on variance.
#' @examples
#' # Example usage of subTopVar function
#' data(iris)
#' dt <- as.data.table(iris)
#' topVarData <- subTopVar(dt, "Species", "Sepal.Length", "Sepal.Width",
#'                         topNperc = 20, complete.rows = TRUE,
#'                         returnLongFormat = FALSE)
#' head(topVarData)
#'
#' @export
subTopVar <- function(dt, row.p, col.p, val, topNperc,
                     min500 = TRUE, topN = NULL,
                     complete.rows, sep.col = "&.&.&",
                     returnLongFormat = FALSE) {
  # priority given to topN
  if (!is.null(topN)) topNperc <- NULL

  datw <- MSqb2:::.l2w(
    dl = dt[, c(..row.p, ..col.p, ..val)],
    col.p = c(col.p),
    row.p = row.p,
    val = val,
    sp = sep.col,
    # fun.agg = median,
    asMat = TRUE,
    Mat.rownm = row.p
  )

  # no missingness
  if (complete.rows) datw <- na.omit(datw)

  # calculate variance
  var.ordered <- apply(datw, 1, function(x) var(x, na.rm = TRUE)) %>%
    sort(., decreasing = TRUE)

  # take top "n" or top "Percent"
  if (!is.null(topNperc)) {
    # take top n% most variant, if this is less than 200 take 200, if data is smaller take everything!
    topN <- round(nrow(datw) * topNperc / 100) %>%
      max(., 200) %>%
      min(., nrow(datw))
    if (min500) {
      topN <- ifelse(topN < 500, 500, topN)
    }
  }
  var.ordered <- names(var.ordered[1:topN])

  # subset topN
  datw <- datw[which(rownames(datw) %in% var.ordered), ]

  MSqb2:::.logg(level = INFO,
               glue("Top {topN} variables with highest variance were used for PCA."))
  
  if (returnLongFormat) {
    return(dt[get(row.p) %in% var.ordered])
  } else {
    return(datw)
  }
}







