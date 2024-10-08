% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msqb_config.R
\name{msqb_config}
\alias{msqb_config}
\title{Configure MSqb2 Analysis Parameters}
\usage{
msqb_config(
  build.para,
  config.file = NULL,
  interactive = TRUE,
  study.variable = "protein",
  filter.by.quaninfo = "auto",
  isolation.interference.cutoff = 75,
  collapse_psm_method = "ionscore",
  filter.singleshot.proteins = "ByPeptide",
  min.intensity = 0.01,
  normalization.method = "vsn",
  normalize.bySubset = NULL,
  summarization.method = "tmp",
  medianpolish.method = "column.eff",
  medianpolish.byPool = FALSE,
  fraction.collapse.method = "maxset",
  batch.corr.method = "none",
  na.imputation.method = "none",
  filter.metadata = NULL,
  feature.annotation.source = "psm.file",
  study.organism = NULL,
  pairwise.contrasts = TRUE,
  manual.contrasts = NULL,
  choose.contrasts = FALSE,
  pairwise.denominator = NULL,
  protein.column = "Master Protein Accessions",
  peptide.seq.column = "Sequence",
  peptide.annot.seq.column = "Annotated Sequence",
  charge.column = "Charge",
  FileID.column = "File ID",
  quaninfo.column = "Quan Info",
  isolation.interference.column = "Isolation Interference [\%]",
  ion.score.column = "Ions Score",
  filename.column = "Spectrum File",
  abundance.columns = "Abundance",
  desc.column = "Master Protein Descriptions",
  contaminant.column = "Contaminant",
  marked.as.column = "Marked as",
  subsymDT = data.table(symout = c("\\\\~", "\\\\+", "\\\\:", "\\\\|", "\\\\-", "\\\\(",
    "\\\\)", "\\\\/", " "), symin = c("_", ".p", ".", ".", "_", ".", ".", ".", ".")),
  ...
)
}
\arguments{
\item{build.para}{A list or object containing the build parameters for the MSqb2 analysis pipeline. This list is generated by the function \code{MSqb_build()}.}

\item{config.file}{A character string specifying the path to a configuration file to load. If \code{NULL}, no configuration file is loaded. Default is \code{NULL}.}

\item{interactive}{Logical. If \code{TRUE}, the function may prompt the user for input during execution. Default is \code{TRUE}.}

\item{study.variable}{A character string specifying the main study variable, either \code{"protein"} or \code{"peptide"}. Default is \code{"protein"}.}

\item{filter.by.quaninfo}{A character string specifying how to filter data based on quantification information. Default is \code{"auto"}.}

\item{isolation.interference.cutoff}{Numeric. The cutoff value for isolation interference, used to filter data. Default is \code{75}.}

\item{collapse_psm_method}{A character string specifying the method for handling multiple features. Default is \code{"ionscore"}.}

\item{filter.singleshot.proteins}{A character string specifying how to filter single-shot proteins. Options are \code{"ByPeptide"} or \code{"ByProtein"}. Default is \code{"ByPeptide"}.}

\item{min.intensity}{Numeric. The minimum intensity threshold for filtering data. Default is \code{0.01}.}

\item{normalization.method}{A character string specifying the normalization method. Options include \code{"vsn"}, \code{"quantile"}, \code{"median"}, etc. Default is \code{"vsn"}.}

\item{normalize.bySubset}{A character string referring to a column name or \code{NULL}. If specified, normalization will be performed by the subset defined by the labels in the specified data column. Default is \code{NULL}.}

\item{summarization.method}{A character string specifying the method for summarizing data. Options include \code{"tmp"}, \code{"median"}, \code{"mean"} and \code{"sum"} Default is \code{"tmp"}.}

\item{medianpolish.method}{A character string specifying the method for median polishing. Options include \code{"column.eff"} and \code{"residual"}. Default is \code{"column.eff"}.}

\item{medianpolish.byPool}{Logical. If \code{TRUE}, median polishing is performed by pool (in TMT method). Default is \code{FALSE}.}

\item{fraction.collapse.method}{A character string specifying the method for collapsing fractions. Options are \code{"MAXset"}, \code{"max"}, \code{"median"}, \code{"sum"}, \code{"single"}, \code{"multi"}. Default is \code{"MAXset"}.}

\item{batch.corr.method}{A character string specifying the method for batch correction. Options include \code{"limma"}, \code{"svn"} and \code{"none"}. Default is \code{"none"}.}

\item{na.imputation.method}{A character string specifying the method for imputing missing data. Options include all methods from packages \code{msImpute} and \code{imputeQRILC.} Default is \code{"none"}.}

\item{filter.metadata}{A character string or \code{NULL}. If specified, metadata will be filtered accordingly. For details see the help of \code{filter_metadata()} function. Default is \code{NULL}.}

\item{feature.annotation.source}{A character string specifying the source of feature annotation. Options include \code{"psm.file"}, \code{"uniprot"}, \code{"uniprot.file"} and \code{"none"}. Default is \code{"psm.file"}.}

\item{study.organism}{A character string specifying the organism studied. Only if \code{feature.annotation.source} is set to \code{"uniprot"}. Default is \code{NULL}.}

\item{pairwise.contrasts}{Logical. If \code{TRUE}, pairwise contrasts are automatically generated. Default is \code{TRUE}.}

\item{manual.contrasts}{A list specifying manual contrasts to be applied. Default is \code{NULL}.}

\item{choose.contrasts}{Logical. If \code{TRUE}, allows the user to choose contrasts interactively. Default is \code{FALSE}.}

\item{pairwise.denominator}{A character string or \code{NULL} specifying the denominator for pairwise contrasts. Default is \code{NULL}.}

\item{protein.column}{A character string specifying the column name for proteins in the input data. Default is \code{"Master Protein Accessions"}.}

\item{peptide.seq.column}{A character string specifying the column name for peptide sequences in the input data. Default is \code{"Sequence"}.}

\item{peptide.annot.seq.column}{A character string specifying the column name for annotated peptide sequences. Default is \code{"Annotated Sequence"}.}

\item{charge.column}{A character string specifying the column name for charge in the input data. Default is \code{"Charge"}.}

\item{FileID.column}{A character string specifying the column name for file IDs in the input data. Default is \code{"File ID"}.}

\item{quaninfo.column}{A character string specifying the column name for quantification information. Default is \code{"Quan Info"}.}

\item{isolation.interference.column}{A character string specifying the column name for isolation interference. Default is \code{"Isolation Interference [\%]"}.}

\item{ion.score.column}{A character string specifying the column name for ion scores. Default is \code{"Ions Score"}.}

\item{filename.column}{A character string specifying the column name for filenames in the input data. Default is \code{"Spectrum File"}.}

\item{abundance.columns}{A character string specifying the column name for abundance values. Default is \code{"Abundance"}.}

\item{desc.column}{A character string specifying the column name for protein descriptions. Default is \code{"Master Protein Descriptions"}.}

\item{contaminant.column}{A character string specifying the column name for contaminants. Default is \code{"Contaminant"}.}

\item{marked.as.column}{A character string specifying the column name for marked entries in the input data. Default is \code{"Marked as"}.}

\item{subsymDT}{A \code{data.table} defining symbolic substitutions, with columns \code{symout} and \code{symin} for output and input symbols respectively.}

\item{...}{Additional arguments passed to other methods or functions.}
}
\value{
A list or object containing the configured parameters for MSqb2 analysis.
}
\description{
This function sets up and configures parameters for an MSqb2 analysis. It allows users to specify various settings related to filtering, normalization, summarization, fraction combination, batch correction, and missing data imputation, among others.
}
\details{
The \code{msqb_config} function provides a comprehensive interface for configuring an MSqb2 analysis. Users can specify various filtering criteria, normalization methods, summarization techniques, and other parameters necessary for processing proteomics data. The function is highly customizable and can be used interactively or with predefined configurations.
}
