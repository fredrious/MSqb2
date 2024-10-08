% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{.find.file}
\alias{.find.file}
\title{Find Files Matching Input Patterns}
\usage{
.find.file(input, output, toData, where)
}
\arguments{
\item{input}{A named list where names are the patterns to search for in file names.}

\item{output}{A named list where names are the expected outputs. The function will append the found file paths to this list.}

\item{toData}{A named list of additional data to be included in the environment after searching.}

\item{where}{A character string specifying the directory where the search should be performed.}
}
\description{
This function searches for files in a specified directory that match patterns provided in the \code{input} parameter.
It logs an error if multiple matches are found for any pattern and updates the environment with the results.
}
\details{
The function searches for files in the \code{where} directory that match each pattern in the \code{input} list.
If multiple files match a single pattern, an error is logged using the \code{.logg} function. The function
appends the file paths to the \code{output} list and updates the environment with \code{output}, \code{input}, and \code{toData}.
}
\examples{
\dontrun{
# Example usage
input_patterns <- list(pattern1 = "file1.csv", pattern2 = "file2.csv")
output_list <- list()
additional_data <- list()
directory <- tempdir()

# Create example files in the directory
file.create(file.path(directory, "file1.csv"))
file.create(file.path(directory, "file2.csv"))

.find.file(input = input_patterns, output = output_list, toData = additional_data, where = directory)

# Check the results in the environment
print(output_list)
}
}
