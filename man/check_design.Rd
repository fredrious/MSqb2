% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_design.R
\name{check_design}
\alias{check_design}
\title{Check Design of Experimental Data}
\usage{
check_design(dat)
}
\arguments{
\item{dat}{A data.table containing experimental design metadata. It should include columns for \code{SampleID}, \code{Condition}, and optionally \code{Block} and \code{TechRep}.}
}
\value{
A character string describing the type of experimental design detected and how the relevant factors are going to be handeled in the statistical model. Possible values include:
\itemize{
\item "fixBlock.fixTechRep" - block: fixed-effect, TechRep: fixed-effect
\item "fixBlock.randomTechRep" - block: fixed-effect, TechRep: random-effect
\item "randomBlock.avgTechRep" - block: random-effect, average over TechRep.
\item "fixBlock" - block: fixed-effect.
\item "randomBlock" - block: random-effect.
\item "fixTech" - TechRep: fixed-effect, no blocking.
\item "randomTech" - TechRep: random-effect, no blocking.
\item "group.comparison" - Simple group comparison design.
}
}
\description{
This function verifies the design of experimental data provided in a data table. It checks for issues such as duplicate column names, multiple condition levels assigned to the same sample, and consistency in block and technical replicate assignments. Based on these checks, it determines the type of experimental design (e.g., fixed block, random block, or group comparison).
}
\examples{
# Example data.table
dat <- data.table(SampleID = c("S1", "S2", "S3", "S4"),
                  Condition = c("C1", "C1", "C2", "C2"),
                  Block = c("B1", "B1", "B2", "B2"),
                  TechRep = c("T1", "T2", "T1", "T2"))
check_design(dat)

}
