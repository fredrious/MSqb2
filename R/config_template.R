config_template <- function() {
## start config (don't touch this!)

  #### --------------- Filters
  study.variable <- "Protein"
  filter.by.quaninfo <- "auto"                       # Character vector. Label(s) in the 'Quan Info' column (if provided) to be used to filter features. If set to 'auto' (default), all features whose 'Quan Info' labels are neither empty nor set to 'unique' will be filtered. If set to 'none', filter will not be applied.
  isolation.interference.cutoff <- 75                # Numeric or NULL. Permitted isolation interference. If numeric, features with an isolation interference larger than 'isolation.interference.cutoff' will be filtered. Default is 50, meaning that features with an isolation intereference larger than 50% will be filtered. If set to NULL or 100, the filter will not be applied.
  multi.features.method <- "IonScore"                       # Character. The arguments sets the strategy how to deal with cases in which there are multiple intensity values available for a given feature per sample and run. Options are 'IonScore' (default), 'mean', 'median', 'max' and 'all'. For more details see the help of function 'MultiPSM'.
  filter.singleshot.proteins <- "ByPeptide"           # Character. Remove single hit proteins. Options are 'ByPeptide' (default) and 'ByFeature'.
  min.intensity <- 0.01                         # Numeric. Intensity threshold below which the intensity values considered to be missing. Default is 0.01.
  reference.channel <- NULL                                # Name of the reference channel in the metadata file. Only relevant for MS method TMT, and if a reference channel exists.
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
