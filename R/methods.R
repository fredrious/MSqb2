
# %%%%%%%%%%%%%%%%%%%%%%%%
#' Import Objects from Parent Environment
#'
#' This function imports all objects from the parent environment of the caller
#' into the current environment. It retrieves the objects from the environment
#' that is two levels up from the current environment and places them into the
#' environment that is one level up.
#'
#' @details The function uses `mget` to get the objects from the parent environment
#' and `list2env` to import them into the current environment. It is useful when
#' you need to bring objects from a higher level environment into the current scope.
#'
#' @return NULL
#' @examples
#' # Example usage:
#' a <- 1
#' b <- 2
#' {
#'   c <- 3
#'   impTopEnv()
#'   # Now `a` and `b` are available in this environment
#'   print(a)  # Should print 1
#'   print(b)  # Should print 2
#' }
#'
#' @export
impTopEnv <- function() {
   list2env(mget(names(parent.frame(n = 2)), envir = parent.frame(n = 2)),
            envir = parent.frame(n = 1))
}


# %%%%%%%%%%%%%%%%%%%%%%%%
#' Export Objects to Parent Environment
#'
#' This function exports all objects from the current environment to the parent
#' environment. It retrieves the objects from the current environment and places
#' them into the environment that is one level up from the current environment.
#'
#' @details The function uses `mget` to get the objects from the current environment
#' and `list2env` to export them to the parent environment. This is useful when you
#' need to make objects available to a higher-level scope.
#'
#' @return NULL
#' @examples
#' # Example usage:
#' {
#'   a <- 1
#'   b <- 2
#'   # a and b are in the current environment
#'   expTopEnv()
#'   # a and b should now be available in the parent environment
#'   print(parent.frame()$a)  # Should print 1
#'   print(parent.frame()$b)  # Should print 2
#' }
#'
#' @export
expTopEnv <- function() {
   list2env(mget(names(parent.frame(n = 1)), envir = parent.frame(n = 1)),
            envir = parent.frame(n = 2))
}





# %%%%%%%%%%%%%%%%%%%%%%%%
#' Expand List into Parent Environment
#'
#' This function takes a named list of objects and imports these objects into
#' the parent environment of the function call. The objects in the list are
#' retrieved using their names and placed into the environment one level up
#' from the current environment.
#'
#' @param x A named list where each element will be imported into the parent environment.
#'
#' @details The function uses `mget` to retrieve objects from the environment
#' specified by the list `x` and `list2env` to export them to the parent environment.
#' This is useful when you want to make multiple objects from a list available
#' in the environment where the function was called.
#'
#' @return NULL
#' @examples
#' # Example usage:
#' a <- 1
#' b <- 2
#' my_list <- list(a = a, b = b)
#' {
#'   c <- 3
#'   expandList(my_list)
#'   # Now `a` and `b` are available in this environment
#'   print(a)  # Should print 1
#'   print(b)  # Should print 2
#' }
#'
#' @export
expandList <- function(x) {
  list2env(mget(names(x), envir = as.environment(x)),
           envir = parent.frame())
}






# %%%%%%%%%%%%%%%%%%%%%%%%
#' Custom Assignment Operator for Parent Environment
#'
#' This custom operator `%<-1%` assigns a value to a variable in the parent
#' environment of the function where it is called. It is a shorthand for
#' assigning a value to an object one level up in the environment hierarchy.
#'
#' @param x A character string specifying the name of the variable to which the
#' value should be assigned.
#' @param y The value to be assigned to the variable.
#'
#' @details The operator uses `assign` to place the value `y` into a variable named `x`
#' in the parent environment of the caller, specifically in the environment that is
#' two levels up from the current environment.
#'
#' @return The value `y`, invisibly.
#' @examples
#' # Example usage:
#' a <- 1
#' {
#'   b %<-1% 2
#'   print(b)  # Should print 2, as `b` is assigned in the parent environment
#' }
#'
#' @export
`%<-1%` <- function(x, y) {
  assign(x, y, pos = parent.frame(n = 2))
}



# %%%%%%%%%%%%%%%%%%%%%%%%
#' Custom Assignment Operator for Global Environment
#'
#' This custom operator `%<-g%` assigns a value to a variable in the global
#' environment, regardless of the current environment where it is called.
#' It provides a convenient way to assign values to variables in the global
#' environment directly from within functions or other environments.
#'
#' @param x A character string specifying the name of the variable to which the
#' value should be assigned.
#' @param y The value to be assigned to the variable.
#'
#' @details The operator uses `assign` to place the value `y` into a variable named `x`
#' in the global environment. This can be useful for modifying global variables
#' from within functions or other non-global environments.
#'
#' @return The value `y`, invisibly.
#' @examples
#' # Example usage:
#' my_global_var <- NULL
#' {
#'   "my_global_var" %<-g% 10
#'   print(my_global_var)  # Should print 10, as `my_global_var` is assigned in the global environment
#' }
#'
#' @export
`%<-g%` <- function(x, y) {
  assign(x, y, pos = .GlobalEnv)
}




# %%%%%%%%%%%%%%%%%%%%%%%%
#' Validate Build Arguments
#'
#' This function performs validation checks on several build arguments. It ensures
#' that the arguments conform to expected types, lengths, and values. The function
#' utilizes various assertion functions to validate the arguments and log assertions
#' using the `MSqb2` package's logging functions.
#'
#' @details The function imports objects from the parent environment using `impTopEnv`
#' and then performs validation on the following arguments:
#' \itemize{
#'   \item `measurements.file` - Must be a character vector.
#'   \item `metadata.file` - Must be a character vector of length 1.
#'   \item `measurements.file.sheet` - Either a character vector of length 1 or an integer, or NULL.
#'   \item `metadata.file.sheet` - Either a character vector of length 1 or an integer, or NULL.
#'   \item `ms.software` - Must be one of the specified choices ("MQ", "PD").
#'   \item `analysis.name` - Character vector of length 1 or NULL.
#'   \item `ms.software` - Character vector of length 1, not NULL.
#'   \item `prefix` - Character vector of length 1 or NULL.
#'   \item `suffix` - Character vector of length 1 or NULL.
#'   \item `add.date.tag` - Logical value.
#' }
#'
#' @return NULL
#' @examples
#' # Example usage:
#' # Assuming these variables are defined and imported from the parent environment:
#' measurements.file <- "data.csv"
#' metadata.file <- "metadata.csv"
#' measurements.file.sheet <- "Sheet1"
#' metadata.file.sheet <- "Sheet1"
#' ms.software <- "MQ"
#' analysis.name <- "Analysis1"
#' prefix <- "prefix_"
#' suffix <- "suffix_"
#' add.date.tag <- TRUE
#'
#' .checkBuildArgs()
#'
#' @import checkmate
#' @export
.checkBuildArgs <- function() {

  impTopEnv() # import all objects from 1 env backward (parent.frame)

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
#' Log Messages with Levels and Outputs
#'
#' This function logs messages with different levels to both a file and the console.
#' It uses the `logger` package to handle the formatting and output.
#'
#' @param level The log level (e.g., "info", "warn", "error", "fatal").
#' @param msg The message to log.
#'
#' @return NULL
#' @import logger magrittr
#' @export
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
#' Create Subdirectories with Custom Naming
#'
#' This function creates a set of subdirectories within a specified path.
#' The names of these subdirectories are constructed using optional prefix,
#' suffix, analysis name, and date tag parameters. It also assigns the
#' paths of these subdirectories to global environment variables for later use.
#'
#' @param path A character string specifying the base path where subdirectories will be created.
#' @param sub.dir A character vector of subdirectory names to be created within the base path.
#' @param prefix An optional character string to prefix each subdirectory name. Default is `NULL`.
#' @param analysis.name An optional character string to include in the subdirectory name. Default is `NULL`.
#' @param suffix An optional character string to suffix each subdirectory name. Default is `NULL`.
#' @param add.date.tag An optional character string representing a date tag to be included in the subdirectory name. Default is `NULL`.
#'
#' @details
#' The function constructs each subdirectory name by concatenating the prefix, subdirectory name, analysis name,
#' suffix, and date tag, separated by underscores. If any of these parameters are `NULL`, they are omitted
#' from the name construction. After creating the directories, the function assigns their paths to global
#' environment variables named `<sub_dir>.path` for each subdirectory.
#'
#' @return NULL
#' @examples
#' # Create subdirectories with specific naming components
#' .mk.dir(
#'   path = tempdir(),
#'   sub.dir = c("data", "results"),
#'   prefix = "project",
#'   analysis.name = "analysis1",
#'   suffix = "v1",
#'   add.date.tag = Sys.Date()
#' )
#'
#' # Check the paths assigned to global environment variables
#' print(data.path)
#' print(results.path)
#'
#' @export
.mk.dir <- function(path, sub.dir, prefix = NULL, analysis.name = NULL,
                    suffix = NULL, add.date.tag = NULL) {
  for (i in seq_along(sub.dir)) {
    sdir.p <- file.path(path,
                        paste(c(prefix, sub.dir[i], analysis.name, suffix, add.date.tag),
                              collapse = "_"))
    if (!dir.exists(sdir.p)) {
      dir.create(sdir.p, recursive = TRUE)
    }
    assign(paste0(sub.dir[i], ".path"), sdir.p, envir = .GlobalEnv)
  }
}




# %%%%%%%%%%%%%%%%%%%%%%%%
#' Find Files Matching Input Patterns
#'
#' This function searches for files in a specified directory that match patterns provided in the `input` parameter.
#' It logs an error if multiple matches are found for any pattern and updates the environment with the results.
#'
#' @param input A named list where names are the patterns to search for in file names.
#' @param output A named list where names are the expected outputs. The function will append the found file paths to this list.
#' @param toData A named list of additional data to be included in the environment after searching.
#' @param where A character string specifying the directory where the search should be performed.
#'
#' @details
#' The function searches for files in the `where` directory that match each pattern in the `input` list.
#' If multiple files match a single pattern, an error is logged using the `.logg` function. The function
#' appends the file paths to the `output` list and updates the environment with `output`, `input`, and `toData`.
#'
#' @return NULL
#' @examples
#' # Example usage
#' input_patterns <- list(pattern1 = "file1.csv", pattern2 = "file2.csv")
#' output_list <- list()
#' additional_data <- list()
#' directory <- tempdir()
#'
#' # Create example files in the directory
#' file.create(file.path(directory, "file1.csv"))
#' file.create(file.path(directory, "file2.csv"))
#'
#' .find.file(input = input_patterns, output = output_list, toData = additional_data, where = directory)
#'
#' # Check the results in the environment
#' print(output_list)
#'
#' @import glue

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
#' Check and Validate File Paths
#'
#' This function verifies the existence of files specified in the `input` list. It checks if these files are
#' either full paths that exist, or if they can be found in a specified directory (`dpth`) or the working directory.
#' If files are missing, the function logs an error or warning and optionally copies the files to the `dpth` directory.
#'
#' @param input A named list where names are file names or patterns and values are file paths. If a value is `NULL`,
#'   it will be ignored.
#' @param dpth A character string specifying the directory path to search for files. Default is `Data.path`.
#' @param toData A logical value indicating whether to copy files to the `dpth` directory if they are not already there. Default is `FALSE`.
#' @param interactive A logical value indicating whether the function should prompt the user interactively to copy files. Default is `TRUE`.
#'
#' @details
#' The function first checks if the files in `input` are full paths and exist. If so, it updates the `output` list
#' and sets `toData` to `TRUE`. If not, it searches for these files in the specified `dpth` directory and, if necessary,
#' in the working directory. If files are still not found, it logs a fatal message.
#' If `toData` is `TRUE` and `interactive` is `TRUE`, the user will be prompted to copy missing files to the `dpth` directory.
#'
#' @return A named list of files with their paths. If files were copied, their paths in the `dpth` directory are returned.
#' @examples
#' # Example usage
#' input_files <- list(file1 = "path/to/file1.csv", file2 = "path/to/file2.csv")
#' output_files <- .check.file(input = input_files, dpth = tempdir(), toData = TRUE, interactive = TRUE)
#' print(output_files)
#'
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
#' Read Data from File
#'
#' Reads data from a file, supporting both Excel files (`.xlsx`, `.xls`) and tab-delimited text files.
#' The function converts the data into a `data.table` object and applies a conversion to factors if necessary.
#'
#' @param file A character string specifying the file path to read. The file can be an Excel file (with `.xlsx` or `.xls` extension)
#'   or a tab-delimited text file.
#' @param sheet An integer specifying the sheet number to read from the Excel file. Default is `1`.
#'
#' @details
#' The function checks the file extension to determine whether it should use `readxl::read_excel` for Excel files or
#' `data.table::fread` for text files. After reading the data, the function converts it to a `data.table` and applies
#' `MSqb2::char2fact` to ensure that character columns are converted to factors.
#'
#' @return A `data.table` object containing the data read from the file.
#' @examples
#' # Example usage with a tab-delimited file
#' data <- read.file("path/to/datafile.txt")
#'
#' # Example usage with an Excel file
#' data <- read.file("path/to/datafile.xlsx", sheet = 2)
#'
#' @import readxl
#' @export
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
#' Import Configuration Parameters from a File
#'
#' Imports configuration parameters from a specified file or generates a new configuration file if none exists.
#' This function is typically called from `msqb_config` and supports both workflow and visualization configurations.
#'
#' @param conf.fl A character string specifying the path to the configuration file. If `NULL`, a new configuration file is generated.
#' @param Scripts.path A character string specifying the path where the scripts and configuration files are located.
#' @param analysis.name A character string specifying the analysis name to include in the configuration file name. Default is `NULL`.
#' @param add.date.tag A character string specifying a date tag to append to the configuration file name. Default is `NULL`.
#' @param whichConf A character string specifying the type of configuration ("workflow" or "viz").
#' @param conf.args A list of additional configuration arguments that can override parameters in the configuration file.
#'
#' @details
#' If the specified configuration file does not exist or is `NULL`, the function creates a new configuration file
#' with default parameters in the specified `Scripts.path`. The function reads the configuration file and evaluates
#' its content in the parent environment.
#'
#' @return A character string specifying the path to the imported or generated configuration file.
#' @import glue
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
#' Log and Assert Conditions
#'
#' Asserts conditions and logs error messages if any assertions fail. This function captures and processes
#' assertion expressions, evaluates them, and logs errors using the `.logg` function if the assertions are not met.
#'
#' @param xpr An expression containing the assertion(s) to evaluate. The expression should be a call to a checkmate
#'   assertion function (e.g., `assertCharacter`, `assertLogical`) without the `add` parameter.
#'
#' @details
#' The function captures the provided expression, modifies it to include the `add` parameter for collecting assertion errors,
#' and then evaluates the modified expression. If any assertion fails, the error messages are collected and logged using the
#' `.logg` function at the error level.
#'
#' @return The function does not return a value. It stops execution and logs an error if any assertion fails.
.loggAssert <- function(xpr) {
   impTopEnv()
   .collAssert = checkmate::makeAssertCollection()
   xpr <- deparse(substitute(xpr)) %>% gsub(" ", "", .) %>% paste(., collapse = "")
   xpr <- paste0(substr(xpr, 1, nchar(xpr)-1), ",add=.collAssert)")
   eval(parse(text = xpr))
   if (!.collAssert$isEmpty()) {
      MSqb2:::.logg(error, skip_formatter(.collAssert$getMessages()))
   }
}




# %%%%%%%%%%%%%%%%%%%%%%%%
#' Convert Factor Columns to Character Columns in a Data Table
#'
#' Converts all factor columns in a data.table to character columns.
#'
#' @param dt A `data.table` object in which factor columns need to be converted to character columns.
#'
#' @details
#' The function identifies all columns in the `data.table` that are of type factor and converts them to character columns.
#' If no factor columns are present, the original `data.table` is returned unchanged.
#'
#' @return A `data.table` object with factor columns converted to character columns.
#'
#' @examples
#' library(data.table)
#' dt <- data.table(a = factor(c("x", "y", "z")), b = 1:3)
#' dt <- fact2char(dt)
#' str(dt) # 'a' column is now character
#'
#' @export
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
#' @param dt Input data of class data.table or data.frame.
#' @return The input data \code{dt} with columns of class character converted to factor.
#' @seealso \code{\link{int2fact}}
#' @examples
#' DT <- data.table(x = c(1., 2., 3.), y = c(4L, 5L, 6L), z = c("a", "b", "NA"))
#' # DT <- data.frame(x = c(1.,2.), y = c(1L,2L), z = c("1","2"))
#' sapply(DT, class)
#' MSqb2::char2fact(DT)
#' sapply(DT, class)
#' @export
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
#' @export
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
#' Match Model Formula with Metadata
#'
#' Checks whether all parameters used in a model formula match the column names in a given metadata design table. If any parameters do not match, the user is prompted to either stop the process and correct the issue or proceed by removing the missing parameters from the model formula.
#'
#' @param dsgn A data frame or data table containing the metadata design information.
#' @param frm A character string representing the model formula to be matched with the metadata design.
#'
#' @details
#' The function extracts all variables used in the model formula and compares them with the column names in the metadata design. If any variables are missing, a message is logged, and the user is prompted to choose whether to stop the process or remove the missing variables from the formula and continue. This function is useful in cases where the metadata design does not fully align with the model formula.
#'
#' @return The updated model formula as a character string, with any missing parameters removed if the user chooses to proceed.
#'
#' @examples
#' dsgn <- data.frame(sampleID = 1:5, group = factor(c("A", "B", "A", "B", "A")))
#' frm <- "~ group + batch"
#' .match_ModelFormula_metadata(dsgn, frm)
#'
#' @import glue
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
#' Read or Generate Configuration Parameters
#'
#' Reads configuration parameters from a provided file or generates them based on build parameters. Handles cases where both configuration parameters and a configuration file are provided, prompting the user to choose which to use. The function can work for both workflow and visualization configuration types.
#'
#' @param build.para A list containing build parameters. This is generated by the `msqb_build` function.
#' @param config.para A list containing pre-existing configuration parameters (optional).
#' @param config.file A character string specifying the path to a configuration file (optional).
#' @param config.type A character string indicating the type of configuration: either `"wf"` for workflow or `"viz"` for visualization.
#' @param ... Additional arguments passed to the configuration functions (`msqb_config` or `msqb_config_viz`).
#'
#' @details
#' This function is designed to handle the reading and generation of configuration parameters based on the build parameters and the specified configuration type. If both `config.file` and `config.para` are provided, the user is prompted to choose which to use. If neither is provided, the function attempts to generate the configuration parameters using the specified configuration function.
#'
#' The function expands the configuration list and exports it to the top environment, making it available for further processing in the workflow.
#'
#' @return The function does not return a value directly. Instead, it modifies the environment by expanding and exporting the configuration parameters.
#'
#' @examples
#' # Assuming `build.para` is already defined
#' readConfigPara(build.para = build.para, config.file = "path/to/config.R", config.type = "wf")
#'
#' @export
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
   expTopEnv()
}




# %%%%%%%%%%%%%%%%%%%%%%%%
#' Convert Wide Data to Long Format
#'
#' Converts a wide-format data table or matrix into long format, where specific columns represent row and column properties. The function also allows splitting identifiers based on a separator and reshapes the data accordingly.
#'
#' @param dw A `data.table` or `matrix` in wide format to be converted.
#' @param val A character string specifying the name of the value column in the output.
#' @param col.p A character vector specifying the column names in the output.
#' @param row.p A character vector specifying the column names in the output.
#' @param sp A character string used as a separator to split the identifiers for row and column properties. Defaults to `"&.&.&"`.
#'
#' @details
#' This function reshapes wide-format data (such as a matrix) into long format, splitting row and column properties based on a specified separator. It is useful for preparing data for downstream analysis, where each unique combination of row and column properties is represented as a separate row in the output.
#'
#' The function handles both matrices and `data.table` objects, allowing for flexibility in input types. Row and column properties are split based on the specified separator and are expanded into new columns.
#'
#' @return A `data.table` in long format with the specified row and column properties, and a value column.
#'
#' @examples
#' # Example with a matrix input
#' matrix_data <- matrix(1:9, nrow = 3, dimnames = list(c("A&.&.&1", "B&.&.&2", "C&.&.&3"), c("X&.&.&1", "Y&.&.&2", "Z&.&.&3")))
#' .w2l(matrix_data, val = "Value", col.p = "Group", row.p = "Category")
#'
#' @import data.table
#' @export
.w2l <- function(dw,
                 val,
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
#' Convert Long Data to Wide Format
#'
#' Converts a long-format `data.table` into a wide-format data table or matrix, where specific columns represent row and column properties. The function also allows combining multiple identifiers using a specified separator.
#'
#' @param dl A `data.table` in long format to be converted.
#' @param col.p A character vector specifying the column properties to be used in the wide format.
#' @param row.p A character vector specifying the row properties to be used in the wide format.
#' @param val A character string specifying the name of the value column in the input.
#' @param sp A character string used as a separator to combine multiple identifiers for row and column properties. Defaults to `"&.&.&"`.
#' @param asMat Logical; if `TRUE`, returns the result as a matrix. Defaults to `FALSE`.
#' @param Mat.rownm A character vector specifying the columns to be used as row names in the matrix output.
#' @param ... Additional arguments passed to `dcast.data.table`, such as `fun.aggregate`.
#'
#' @details
#' This function reshapes long-format data into wide format, combining row and column properties using the specified separator. It is useful for preparing data for downstream analysis, where each unique combination of row and column properties is represented as a separate column.
#'
#' The function handles both `data.table` objects and allows the output to be returned as either a `data.table` or a matrix. If a matrix is returned, optional row names can be specified.
#'
#' @return A wide-format `data.table` or matrix.
#'
#' @examples
#' # Example with a data.table input
#' library(data.table)
#' long_data <- data.table(Category = rep(c("A", "B", "C"), each = 3),
#'                         Group = rep(c("X", "Y", "Z"), 3),
#'                         Abundance = 1:9)
#' wide_data <- .l2w(long_data, col.p = "Group", row.p = "Category")
#'
#' # Example returning a matrix
#' wide_matrix <- .l2w(long_data, col.p = "Group", row.p = "Category", asMat = TRUE, Mat.rownm = "Category")
#'
#' @import data.table
#' @export
.l2w <- function(dl,
                 col.p,
                 row.p,
                 val,
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
#' Clean and Standardize Data Column Names and Values
#'
#' This function cleans and standardizes the column names and optionally the data within a `data.table` or vector by replacing or removing symbols that may interfere with design formulas. It is particularly useful in preprocessing data for statistical analysis where certain characters in column names or data values can cause issues.
#'
#' @param dt A `data.table` or vector containing the data to be cleaned.
#' @param complete.check Logical; if `TRUE`, the function will also clean the content within the columns of the `data.table`. Defaults to `FALSE`.
#' @param subsymDT A `data.table` specifying the symbols to be replaced or removed. It should contain two columns: `"symout"` (the symbol to be replaced) and `"symin"` (the replacement). Defaults to a table with common problematic symbols.
#'
#' @details
#' The function ensures that certain symbols commonly used in design formulas (e.g., `~`, `+`, `:`, `|`, `-`, `(`, `)`, `/`, and spaces) are either removed or replaced in column names and, if `complete.check = TRUE`, within the data itself. This avoids potential conflicts when working with formulas in statistical models.
#'
#' The default `subsymDT` table maps symbols to replacements:
#' \itemize{
#'   \item `"~"` -> `""`
#'   \item `"+"` -> `"."`
#'   \item `":"` -> `""`
#'   \item `"|"` -> `""`
#'   \item `"-"` -> `"_"`
#'   \item `"("` -> `"._"`
#'   \item `")"` -> `"_."`
#'   \item `"/"` -> `"."`
#'   \item `" "` -> `""`
#' }
#'
#' @return A list with two elements:
#' \itemize{
#'   \item `"cln.nm.dt"`: The cleaned `data.table` or vector.
#'   \item `"ref.name.tbl"`: A reference table showing the original and cleaned column names (if `dt` is a `data.table`).
#' }
#'
#' @examples
#' library(data.table)
#' dt <- data.table(`Column (1)` = c(1, 2), `Another+Column` = c(3, 4))
#' result <- cleanDataPara(dt)
#' cleaned_data <- result$cln.nm.dt
#' reference_table <- result$ref.name.tbl
#'
#' @import limma
#' @export

cleanDataPara <- function(dt, complete.check = FALSE, # export2env = FALSE
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



  if (is.data.table(dt)) {
    nmDT <- data.table(original.names = names(dt)) # old names with punctuations
    for (is in seq(nrow(subsymDT))) {
      names(dt) <- gsub(subsymDT$symout[is], subsymDT$symin[is], names(dt))
    }
    names(dt) <- make.names(names(dt))
    names(dt) <- sapply(names(dt), function(x) {
      limma::strsplit2(x, "\\.") %>%
        .[nchar(.) > 0] %>%
        paste(., collapse = ".")
    })
    nmDT[, new.names := names(dt)] # new names with punctuations

    if (complete.check) {
      for (is in seq(nrow(subsymDT))) {
        for (para in nmDT$new.names) {
          dt[, (para) := lapply(.SD, function(x) gsub(subsymDT$symout[is], subsymDT$symin[is], x, perl = TRUE)), .SDcols = para]
        }
      }
    }

  }

  if (is.vector(dt)) {
    for (is in seq(nrow(subsymDT))) {
      dt <- gsub(subsymDT$symout[is], subsymDT$symin[is], dt)
      nmDT <- NULL
    }
  }

  dout <- list("cln.nm.dt" = dt, "ref.name.tbl" = nmDT)

  return(dout)
}




# %%%%%%%%%%%%%%%%%%%%%%%%
#' Order a `data.table` by Alphanumeric Column
#'
#' This function orders a `data.table` based on alphanumeric sorting of a specified column. The sorting is done by extracting numeric components from the column values while preserving the overall alphanumeric order. It is particularly useful for cases where numeric sequences are embedded within text.
#'
#' @param dt A `data.table` containing the data to be ordered.
#' @param x A string specifying the column name to be ordered. The column can be either character or factor type.
#'
#' @details
#' The function extracts numeric values from the specified column and orders the `data.table` based on these numbers. It handles mixed alphanumeric values by focusing on the numeric part while maintaining the overall alphanumeric order. If the column is a factor, the levels are also updated to reflect the new order.
#'
#' The underlying idea is based on a solution from [Stack Overflow](https://stackoverflow.com/questions/49084625/sort-strings-by-numbers-inside-them).
#'
#' @return The input `data.table` ordered by the specified column.
#'
#' @examples
#' library(data.table)
#' dt <- data.table(id = c("item1", "item10", "item2", "item20"))
#' ordered_dt <- order.num(dt, "id")
#' print(ordered_dt)
#'
#' @export
order.num <- function(dt, x) {
  xo <- dt[, ..x] %>% unlist()
  fct <- ifelse(is.factor(xo), TRUE, FALSE)
  if (fct) xo %>% unique() %>% as.character()
  ix <- xo %>% gsub("[^[:digit:]]", "", .) %>% as.numeric() %>% order() %>% xo[.]
  dt <- dt[order(match(get(x), ix))]
  if (fct) dt[, (x) := factor(get(x), levels = ix)]
  return(dt)
}








