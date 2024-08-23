#' Generate a Configuration Template for MSqb2 Analysis
#'
#' This function generates a configuration template for MSqb2 analysis, providing a set of predefined parameters that can be customized for specific mass spectrometry data processing tasks.
#'
#' @return The function does not return any value but sets up a configuration template with predefined parameters. This template includes various settings for filtering, normalization, summarization, fraction combination, batch effect removal, missing data imputation, annotation, contrast definition, and reformatting of input data.
#'
#' @details The `config_template` function initializes a set of configuration parameters for MSqb2 analysis, which can be modified by the user as needed. The parameters are grouped into several categories:
#'
#' \itemize{
#'   \item **Filters**: Settings for filtering data based on quantification information, isolation interference, intensity thresholds, and more.
#'   \item **Normalization**: Parameters for normalizing the data, including method selection and subset normalization.
#'   \item **Summarization**: Options for summarizing data, including the method for median polishing.
#'   \item **Fraction Combination**: Method for combining data across fractions.
#'   \item **Batch Effect Removal**: Method for removing batch effects from the data.
#'   \item **Missing Data Imputation**: Parameters for handling missing data.
#'   \item **Metadata Filtering**: Options for filtering the dataset based on metadata.
#'   \item **Annotation**: Source of feature annotation and the organism under study.
#'   \item **Contrasts**: Parameters for defining contrasts in differential expression analysis, including pairwise contrasts and complex contrast definitions.
#'   \item **Reformatting**: Options for reformatting input data, particularly for data generated from Proteome Discoverer.
#' }
#'
#' The template is intended to be customized by users for their specific analysis needs, ensuring that all relevant parameters are set before running the MSqb2 pipeline.
#'
#' @examples
#' \dontrun{
#'   # Generate a configuration template
#'   config_template()
#'
#'   # Users can modify the parameters as needed before running the analysis
#' }
#'
#' @export
config_template <- function() {
## start config (don't touch this!)

  #### --------------- Filters
  study.variable <- "Protein"
  filter.by.quaninfo <- "auto"                       # Character vector. Label(s) in the 'Quan Info' column (if provided) to be used to filter features. If set to 'auto' (default), all features whose 'Quan Info' labels are neither empty nor set to 'unique' will be filtered. If set to 'none', filter will not be applied.
  isolation.interference.cutoff <- 75                # Numeric or NULL. Permitted isolation interference. If numeric, features with an isolation interference larger than 'isolation.interference.cutoff' will be filtered. Default is 50, meaning that features with an isolation intereference larger than 50% will be filtered. If set to NULL or 100, the filter will not be applied.
  collapse_psm_method <- "IonScore"                       # Character. The arguments sets the strategy how to deal with cases in which there are multiple intensity values available for a given feature per sample and run. Options are 'IonScore' (default), 'mean', 'median', 'max' and 'all'. For more details see the help of function 'collapse_psm'.
  filter.singleshot.proteins <- "ByPeptide"           # Character. Remove single hit proteins. Options are 'ByPeptide' (default) and 'ByFeature'.
  min.intensity <- 0.01                         # Numeric. Intensity threshold below which the intensity values considered to be missing. Default is 0.01.
  #### --------------- Normalization
  normalization.method = "vsn"
  normalize.bySubset = NULL
  #### --------------- Summarization
  summarization.method = "tmp"
  medianpolish.method = "column.eff"
  medianpolish.byPool = FALSE
  #### ---------------Ftaction combination
  fraction.collapse.method = "maxset"
  #### --------------- remove batch effect
  batch.corr.method = "none"
  #### --------------- impute missing data
  na.imputation.method = "none"
  #### --------------- Filter metadata
  ## (this will consequently filter the dataset)
  filter.metadata = NULL
  #### --------------- Annotation
  feature.annotation.source = "psm.file"
  study.organism = NULL
  #### --------------- defining contrasts
  pairwise.contrasts = TRUE
  manual.contrasts = NULL
  make.complex.contrasts = FALSE
  limma.blocking.parameter = NULL
  limma.trend = TRUE # either of TRUE, FALSE (from limma) or "PSM.count" or "PSM.mean"
  limma.robust = TRUE
  choose.contrasts = FALSE
  pairwise.denominator <- "WT"
  #### ---------------reformat
  ## if input data is PSM-level excel file from 'Proteome Discoverer'
  Sample.col <- "Sample"                            # Name of column containing sample names in metadata.
  protein.column <- "Master Protein Accessions"        # Name of "Master Protein Accessions" column in excel file.
  peptide.seq.column <- "Sequence"               # Name of "Sequence" column in excel file.
  peptide.annot.seq.column <- "Annotated Sequence"               # Name of "Annotated Sequence" column in excel file.
  charge.column <- "Charge"                            # Name of "Charge" column in excel file.
  FileID.column <- "File ID"                            # Name of "File ID" column in excel file.
  quaninfo.column <- "Quan Info"                       # Name of "Quan Info" column in excel file.
  isolation.interference.column <- "Isolation Interference [%]"      # Name of "Isolation Interference" column in excel file.
  ion.score.column <- "Ions Score"                         # Name of "Ions Score" column in excel file.
  filename.column <- "Spectrum File"                   # Name of "Spectrum File" column in excel file.
  abundance.columns <- "Abundance"                      # The prefix term of all columns containing PSM intensity values. For example, if the MS method is TMT, the names of abundance columns usually have the following format: Abundance: 129C, Abundance: 130N, ... . The common term between all channels is then 'Abundance'. This common name should be passed to the argument 'abundance.columns'.
  desc.column <- "Master Protein Descriptions"
  contaminant.column <- "Contaminant"
  marked.as.column <-  "Marked as"



  #### extra parameters ----

## end config (don't touch this!)
}
