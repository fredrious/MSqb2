#' Check and Process Input File for MS Analysis
#'
#' This function checks the input file or data table for MS analysis, determines the data level (PSM, peptide, or protein), identifies relevant columns, and processes the data to extract useful metadata.
#'
#' @param x A `data.table` or a character string representing the path to an input file containing MS data. If a file path is provided, the file will be read into a `data.table`.
#' @param sheet An integer specifying the sheet number to read from if the input file is an Excel file. Default is `1`.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item `"MSmethod"`: The identified mass spectrometry method (`"TMT"` or `"LFQ"`).
#'   \item `"data.level"`: The data level of the input data (`"PSM-level"`, `"Peptide-level"`, or `"Protein-level"`).
#'   \item `"Intensity.col"`: The identified columns corresponding to intensity values.
#'   \item `"Sample.dt"`: A `data.table` containing metadata extracted from the input data, such as file IDs, channels, and conditions.
#'   \item `"parList"`: A list of identified relevant column names for further processing.
#' }
#'
#' @details The `check_input_file` function performs several key tasks:
#' \itemize{
#'   \item Determines whether the data is at the PSM, peptide, or protein level based on the presence of specific columns.
#'   \item Identifies the mass spectrometry method (TMT or LFQ) by examining the intensity columns.
#'   \item Extracts relevant metadata, including sample IDs, conditions, and channels, and organizes it into a structured `data.table`.
#'   \item Prepares a list of relevant column names for use in downstream analyses.
#' }
#'
#' The function supports both direct data input as a `data.table` and reading data from an input file. It is designed to handle data from common MS analysis software (PD and MQ) outputs.
#'
#' @importFrom magrittr %>% %<>%
#'
#' @examples
#' \dontrun{
#'   # Check an input file for MS analysis
#'   input_file <- "path/to/ms_data.xlsx"
#'   result <- check_input_file(input_file, sheet = 1)
#'   print(result)
#' }
#'
#' @export
check_input_file <- function(x, sheet = 1) {

  if (is.character(x)) {
    suppressWarnings( {x <- read.file(file = x, sheet = sheet) %>% unique()} )
  }
  columns <- names(x)


  # data level: psm, peptide or protein
  if (any(c("Intensity", "Ions Score", "XCorr", "Charge") %in% columns)) {

    dt.level <- "PSM-level"

  } else {

    if ("Sequence" %in% columns) {
      if (any(c("Positions in Master Proteins", "Annotated Sequence") %in% columns)) {
        dt.level <- "Peptide-level"
      } else {
        xd <- x[, .(nSeq = uniqueN(Sequence)), by = c("Accession")][nSeq > 1]
        dt.level <- ifelse(nrow(xd) > 1, "Peptide-level", "Protein-level")
      }
    }

    if ("Accession" %in% columns & any(columns %like% "Coverage")) {
      dt.level <- "Protein-level"
    }

  }



  # subset columns to only needed columns
  #######################################
  if (dt.level == "PSM-level") {
    parList <- c(
      protein.column = "Master Protein Accessions",
      peptide.seq.column = "Sequence",
      peptide.annot.seq.column = "Annotated Sequence",
      charge.column = "Charge",
      FileID.column = "File ID",
      ion.score.column = intersect(c("Ions Score", "XCorr"), columns)[1], # forces "Ions Score" if both exists
      quaninfo.column = "Quan Info",
      isolation.interference.column = "Isolation Interference [%]",
      filename.column = "Spectrum File",
      desc.column = "Master Protein Descriptions",
      contaminant.column = "Contaminant")
  }


  if (dt.level == "Peptide-level") { # Peptide-level data
    parList <- c(
      protein.column = "Protein Accessions",
      peptide.seq.column = "Sequence",
      quaninfo.column = "Quan Info",
      contaminant.column = "Contaminant")
  }


  if (dt.level == "Protein-level") { # Protein-level data
    parList <- c(
      protein.column = "Accession",
      contaminant.column = "Contaminant",
      desc.column = "Description",
      gene.column = "Gene Symbol")
  }

  ## 'Marked as' may or may not be in the exported data
  if ("Marked as" %in% columns) parList <- c(parList, c(marked.as.column = "Marked as"))
  #######################################



  # select abundance columns if not psm
  abnd.cols <- columns %>% .[. %like% "Abundance" & !(. %like% "Count|Grouped|Ratio|Normalized|Scaled|Precursor")]


  # check if there are tmt channels in the data, if yes, TMT  = TRUE
  is.tmt.dt <- data.table(cols = abnd.cols) %>%
    .[, cols := gsub("[[:punct:]]", "", cols)] %>%
    .[, tstrsplit(cols, "[[:blank:]]")]

  ix <- apply(is.tmt.dt, 2, \(x) length(intersect(gsub("[[:alpha:]]", "", x), 126:133)) > 3 ) # minimum 5 channels

  if (any(ix)) {
    MSmethod <- ifelse(
      is.tmt.dt[, ix, with = FALSE] %>%
        unlist() %>%
        gsub("[[:digit:]]", "", .) %>%
        .[nchar(.) > 0] %>%
        unique() %>%
        sort() %>%
        identical(., c("C", "N")),
      "TMT", "LFQ")
  } else MSmethod <- "LFQ"

  # generate metadata guide in PSM-level, using file id
  if (dt.level == "PSM-level" & "File ID" %in% columns) mdt <- decompFileID(x, MSmethod)

  # in lfq, psm, use Intensity column
  if (MSmethod == "LFQ" & dt.level == "PSM-level") abnd.cols <- columns[tolower(columns) == "intensity"]



  # if input is protein or peptide level
  if (dt.level %in% c("Peptide-level", "Protein-level")) {

    mdt <- gsub("Abundance|Sample|Control|n/a| ", "", abnd.cols) %>%
      strsplit(., split = ",|\\:") %>%
      do.call(rbind, .) %>%
      as.data.table() %>%
      .[, which(!duplicated( t(.))), with = F] %>%
      .[, sapply(., \(x) all(unique(x) != "")), with = F]


    if (MSmethod == "TMT") {
      cls <- c("FileID", "Channel", "Condition")
    } else {
      cls <- c("FileID", "Condition")
    }

    mdt %<>% setnames(., names(.[, seq_along(cls), with = F]), cls) %>%
      .[, SampleID := apply(.[, -"FileID"], 1, paste, collapse = "_")] %>%
      .[, Sample.cols := abnd.cols] %>%
      .[, c("Sample.cols", "SampleID", cls), with = FALSE]
  }

  cat("\nMSmethod:", MSmethod)
  cat("\nInput data:", dt.level, "\n")
  print(mdt)
  return(invisible( list(MSmethod = MSmethod,
                         data.level = dt.level,
                         Intensity.col = abnd.cols,
                         Sample.dt = mdt,
                         parList = parList) ))
}
