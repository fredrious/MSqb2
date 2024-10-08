% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DE_stat_func.R
\name{DE_stat_func}
\alias{DE_stat_func}
\title{Perform Differential Expression Analysis Using LIMMA}
\usage{
DE_stat_func(
  dat,
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
  pairwise.denominator = c("wt", "wildtype", "wild-type", "cntrl", "control", "ctl",
    "untreated"),
  manual.contrasts = NULL,
  choose.contrasts = FALSE,
  make.complex.contrasts = FALSE,
  blocking.parameter = NULL,
  reverse.contrasts = FALSE,
  showWarnings = FALSE,
  prefix = NULL,
  suffix = NULL,
  rmSingleShotProteins = TRUE,
  add.vs.inContrast = TRUE,
  path = getwd(),
  ...
)
}
\arguments{
\item{dat}{A \code{data.table} containing the MS data. The table should include columns for samples, variables (e.g., proteins),
and values (e.g., abundances).}

\item{sample}{A character string specifying the column name that contains sample identifiers. Default is \code{"SampleID"}.}

\item{variable}{A character string specifying the column name that contains the variables of interest (e.g., \code{"Protein"}). Default is \code{"Protein"}.}

\item{value}{A character string specifying the column name that contains the values to be analyzed (e.g., \code{"Abundance"}). Default is \code{"Abundance"}.}

\item{model.formula}{A character string specifying the model formula to be used in the analysis. Default is \code{"~ 0 + Condition"}.}

\item{sva.obj}{An optional \code{sva} object for batch correction. Default is \code{NULL}.}

\item{techrep}{A character string specifying the column name that contains technical replicates. If provided, these will be used as blocking factors. Default is \code{NULL}.}

\item{limma.trend}{Logical or character. If \code{TRUE}, trend fitting is applied; if a character, it specifies the column name to be used for trend fitting. Default is \code{TRUE}.}

\item{limma.robust}{Logical. If \code{TRUE}, robust fitting is applied in the eBayes step. Default is \code{TRUE}.}

\item{limma.winsor.tail.p}{Logical or numeric. If \code{FALSE}, no winsorization is applied; if numeric, specifies the tail probability for winsorization. Default is \code{FALSE}.}

\item{time.course}{Logical. If \code{TRUE}, a time course analysis is performed using spline regression. Default is \code{FALSE}.}

\item{spline.df}{Integer. Degrees of freedom for spline regression in time course analysis. Default is \code{4}.}

\item{removeTotalMissing}{Logical. If \code{TRUE}, rows with all missing values are removed from the analysis. Default is \code{TRUE}.}

\item{pairwise.contrasts}{Logical. If \code{TRUE}, all pairwise contrasts between conditions are computed. Default is \code{TRUE}.}

\item{pairwise.denominator}{Character vector. Specifies the condition levels to be used as denominators in pairwise contrasts. Default includes \code{"wt"}, \code{"wildtype"}, \code{"control"}, etc.}

\item{manual.contrasts}{A list of manual contrasts to be applied. Default is \code{NULL}.}

\item{choose.contrasts}{Logical. If \code{TRUE}, allows manual selection of contrasts from a list. Default is \code{FALSE}.}

\item{make.complex.contrasts}{Logical. If \code{TRUE}, prompts the user to create complex contrasts. Default is \code{FALSE}.}

\item{blocking.parameter}{A character string specifying the column name for blocking factors, typically technical replicates. Default is \code{NULL}.}

\item{reverse.contrasts}{Logical. If \code{TRUE}, reverses the direction of the contrasts. Default is \code{FALSE}.}

\item{showWarnings}{Logical. If \code{TRUE}, warnings are shown during execution. Default is \code{FALSE}.}

\item{prefix}{A character string to prefix the output filenames. Default is \code{NULL}.}

\item{suffix}{A character string to suffix the output filenames. Default is \code{NULL}.}

\item{rmSingleShotProteins}{Logical. If \code{TRUE}, proteins with only one peptide are removed from the analysis. Default is \code{TRUE}.}

\item{add.vs.inContrast}{Logical. If \code{TRUE}, the contrasts are labeled with "vs." between conditions. Default is \code{TRUE}.}

\item{path}{A character string specifying the output directory for saving results. Default is the current working directory.}

\item{...}{Additional arguments passed to the underlying LIMMA functions.}
}
\value{
A list containing:
\itemize{
\item \code{"fit_eb"}: The fitted model after applying eBayes.
\item \code{"fit_residuals"}: The residuals of the fit.
\item \code{"fit0"}: The initial fitted model before applying contrasts.
}
}
\description{
This function performs differential expression analysis on MS data using limma (Linear Models for Microarray Data) package.
It supports a wide range of experimental designs, including handling of technical replicates, blocking factors, and the application of contrasts
to identify differentially expressed proteins or peptides.
}
\details{
The \code{DE_stat_func} function is designed to handle complex experimental designs typical in mass spectrometry-based proteomics studies.
It supports:

\itemize{
\item \strong{Technical Replicates}: When \code{techrep} is provided, the function accounts for repeated measures by using the \code{duplicateCorrelation} method from LIMMA.
\item \strong{Time Course Analysis}: If \code{time.course} is \code{TRUE}, the function performs spline regression to model the effect of time.
\item \strong{Contrast Definitions}: Users can define contrasts manually or choose from pairwise contrasts. Complex contrasts, such as average or delta-delta contrasts, can also be generated.
}

This function is suitable for both differential expression analysis in standard proteomics experiments and more complex designs involving time courses or technical replicates.
When working with experiments that have potential batch effects, users can incorporate surrogate variable analysis (\code{sva.obj}) to correct for these.
}
