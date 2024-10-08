% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{readConfigPara}
\alias{readConfigPara}
\title{Read or Generate Configuration Parameters}
\usage{
readConfigPara(build.para, config.para, config.file, config.type, ...)
}
\arguments{
\item{build.para}{A list containing build parameters. This is generated by the \code{msqb_build} function.}

\item{config.para}{A list containing pre-existing configuration parameters (optional).}

\item{config.file}{A character string specifying the path to a configuration file (optional).}

\item{config.type}{A character string indicating the type of configuration: either \code{"wf"} for workflow or \code{"viz"} for visualization.}

\item{...}{Additional arguments passed to the configuration functions (\code{msqb_config} or \code{msqb_config_viz}).}
}
\value{
The function does not return a value directly. Instead, it modifies the environment by expanding and exporting the configuration parameters.
}
\description{
Reads configuration parameters from a provided file or generates them based on build parameters. Handles cases where both configuration parameters and a configuration file are provided, prompting the user to choose which to use. The function can work for both workflow and visualization configuration types.
}
\details{
This function is designed to handle the reading and generation of configuration parameters based on the build parameters and the specified configuration type. If both \code{config.file} and \code{config.para} are provided, the user is prompted to choose which to use. If neither is provided, the function attempts to generate the configuration parameters using the specified configuration function.

The function expands the configuration list and exports it to the top environment, making it available for further processing in the workflow.
}
\examples{
\dontrun{
# Assuming `build.para` is already defined
readConfigPara(build.para = build.para, config.file = "path/to/config.R", config.type = "wf")
}

}
