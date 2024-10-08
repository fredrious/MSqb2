% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{.mk.dir}
\alias{.mk.dir}
\title{Create Subdirectories with Custom Naming}
\usage{
.mk.dir(
  path,
  sub.dir,
  prefix = NULL,
  analysis.name = NULL,
  suffix = NULL,
  add.date.tag = NULL
)
}
\arguments{
\item{path}{A character string specifying the base path where subdirectories will be created.}

\item{sub.dir}{A character vector of subdirectory names to be created within the base path.}

\item{prefix}{An optional character string to prefix each subdirectory name. Default is \code{NULL}.}

\item{analysis.name}{An optional character string to include in the subdirectory name. Default is \code{NULL}.}

\item{suffix}{An optional character string to suffix each subdirectory name. Default is \code{NULL}.}

\item{add.date.tag}{An optional character string representing a date tag to be included in the subdirectory name. Default is \code{NULL}.}
}
\description{
This function creates a set of subdirectories within a specified path.
The names of these subdirectories are constructed using optional prefix,
suffix, analysis name, and date tag parameters. It also assigns the
paths of these subdirectories to global environment variables for later use.
}
\details{
The function constructs each subdirectory name by concatenating the prefix, subdirectory name, analysis name,
suffix, and date tag, separated by underscores. If any of these parameters are \code{NULL}, they are omitted
from the name construction. After creating the directories, the function assigns their paths to global
environment variables named \verb{<sub_dir>.path} for each subdirectory.
}
\examples{
\dontrun{
# Create subdirectories with specific naming components
.mk.dir(
  path = tempdir(),
  sub.dir = c("data", "results"),
  prefix = "project",
  analysis.name = "analysis1",
  suffix = "v1",
  add.date.tag = Sys.Date()
)

# Check the paths assigned to global environment variables
print(data.path)
print(results.path)
}

}
