#' Construct the skeleton of the analysis
#'
#' The function checks the input files and construct a folder structure to store output of the workflow.
#' By default, the current directory from which the function is called will be considered as the project
#' path and the basename of the path will be considered as the project name.
#' Three main directories will be created under the project path:
#' 'Data', 'Analysis' and 'Scripts'. The output from each round of the
#' analysis will be saved under 'Analysis' in a separate directory named by the argument 'analysis.name'.
#' If 'analysis.name' is NULL, the output files will be stored directly in the 'Analysis' directory.
#' A log file will be generated for each call of the function 'msqb_build'. The log files have a date
#' and time tag and will be saved in the '_log_files' sub-directory under 'Analysis' directory.
#'
#'
#' @param measurements.file character string or vector, naming the file(s) containing the MS measurements. The input file(s) must either be the PSM-level file(s) from ProteomDiscoverer or evidence file from MaxQuant. If more than one filename provided, the files will be merged to construct a single dataset.
#' @param metadata.file character string, naming the excel file containing the experimental design.
#' @param measurements.file.sheet either a character string (the name of a sheet), or an integer (the position of a sheet). Only relevant if measurements.file is from xlsx family. By default the first sheet will be read.
#' @param metadata.file.sheet either a character string (the name of a sheet), or an integer (the position of the sheet). By default the first sheet will be read.
#' @param ms.software character, naming the software used to generate ms data. Options are 'MQ' for MaxQuant and 'PD' for Proteome Discoverer.
#' @param analysis.name character string. If provided, a directory with the provided name will be created under the Analysis directory. The output of the workflow will then be saved in this directory.
#' @param prefix character string to be added at the beginning of the directory names.
#' @param suffix character string to be added at the end of the directory names.
#' @param add.date.tag logical. If TRUE, the date of the analysis will be added to the Analysis directory. Default is TRUE.
#' @param interactive logical.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom glue glue
#' @return A list containing the paths and input parameters.
#' @export
#'
msqb_build <- function(measurements.file, # either psm file from ProteomDiscoverer or evidence file from MaxQuant
                       metadata.file, # Experimental design
                       measurements.file.sheet = 1, # if measurements.file is from xlsx family
                       metadata.file.sheet = 1, # if metadata.file is from xlsx family
                       ms.software,
                       analysis.name = NULL,
                       prefix = NULL,
                       suffix = NULL,
                       add.date.tag = TRUE,
                       interactive = TRUE) {


   # test input arguements
   .checkBuildArgs()


   # make top dirs
   .mk.dir(path = getwd(), sub.dir = c("Data", "Analysis", "Scripts"))


   # make data tag
   add.date.tag %<>% ifelse(., paste0("_", format(Sys.time(), "%Y%m%d")), "")



   # recreate Analysis.path if necessary!
   if (!unique(c(prefix, analysis.name, suffix, add.date.tag)) %>% is.null()) {
    .mk.dir(Analysis.path, "Analysis", prefix, analysis.name, suffix, add.date.tag)
   }



  # create sub-directories in Analysis.path for saving files!
   .mk.dir(path = Analysis.path,
           sub.dir = c("QCplots", "Volcanos", "PCA", "Tables", "_log_files",
                       "TopSig_Heatmap", "HighVar_Heatmap", "ProfilePlots", "VennDiagram"))


  # analysis logging
  log_file <<- file.path(Analysis.path, "_log_files",
                        paste0("_log_", format(Sys.time(), "%Y%m%d_%Hh%Mm%Ss"), ".log"))
  .logg(info, "Working dir:\n{getwd()}\nProject's name:\n{basename(getwd())}")



  # check input data - copy to 'Data' sub-directory if approved by user!
  input.files <- list("measurements.file" = measurements.file,
                      "metadata.file" = metadata.file) %>%
     .check.file(., dpth = Data.path, toData = interactive, interactive = interactive)
  for (i in seq_along(input.files)) assign(names(input.files[i]), input.files[[i]])



  ## check sheet names when excel file ----
  for (ff in input.files) {
    if ((tolower(strsplit(ff[1], split = "\\.")[[1]][-1]) %in% c("xls", "xlsx"))) {
      if (is.null(paste0(ff, ".sheet"))) {
         .logg(WARN,
            glue( "The sheet name of the input file\n {basename(ff)} \nnot provided. ",
                  "By default the first sheet will be read."))
        assign(paste0(ff, ".sheet"), 1)
      }
    }
  }




  ## check measurements.file ----
  if (!tolower(strsplit(measurements.file[1], split = "\\.")[[1]][-1]) %in% c("xls", "xlsx", "txt")) {
     .logg(fatal,
        glue("The Measurements file must either be from xlsx family or a txt file."))
  }

  rm(input.files, i, ff)
  mget(ls(name = environment()), envir = environment()) %>%
    return(.)
}
