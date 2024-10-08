% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantile_normalisation.R
\name{quantile_normalisation}
\alias{quantile_normalisation}
\title{Quantile Normalization of Data}
\usage{
quantile_normalisation(dat, col.p, row.p, val, log2transform = TRUE)
}
\arguments{
\item{dat}{A \code{data.table} containing the data to be normalized.}

\item{col.p}{Character vector. Column names used for grouping the data (e.g., sample or condition).}

\item{row.p}{Character vector. Row identifiers, typically the feature or protein IDs.}

\item{val}{Character. The column name representing the intensity or abundance values to be normalized.}

\item{log2transform}{Logical. Whether to log2 transform the data before normalization. Default is \code{TRUE}.}
}
\value{
A \code{data.table} with quantile normalized intensity values.
}
\description{
This function performs quantile normalization on mass spectrometry data, with an option to log2 transform the data before normalization.
}
\details{
The function uses the \code{normalize.quantiles} function from the \code{preprocessCore} package to perform quantile normalization. This method makes the distribution of intensity values the same across samples. Log2 transformation is applied before normalization if \code{log2transform} is set to \code{TRUE}.
}
\examples{
\dontrun{
  normalized_data <- quantile_normalisation(dat = my_data, col.p = "SampleID", row.p = "ProteinID", val = "Intensity")
}

}
