% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation_query.R
\name{annotation_query}
\alias{annotation_query}
\title{Annotate Protein Data with Gene Information}
\usage{
annotation_query(
  dat,
  feature.annotation.source,
  uniprot.input = "Protein",
  uniprot.output = "Genes",
  uniprot.extra = NULL,
  uniprot.KEY = "UNIPROTKB",
  uniprot.species.name = NULL,
  uniprot.taxId = NULL,
  n.call = 99,
  uniprot.file = "AnnotationsTable.xlsx",
  shorten.gene.names = TRUE,
  split.char = ";",
  save.annotation = TRUE,
  save.annotation.path = Tables.path,
  Data.path = Data.path
)
}
\arguments{
\item{dat}{A data table containing protein data with the columns specified by \code{uniprot.input}.}

\item{feature.annotation.source}{Character string specifying the source of annotation. Options are \code{"uniprot.file"}, already generated \code{"annotation.file"} from a previous run, or \code{"uniprot"}.}

\item{uniprot.input}{Character string specifying the column in \code{dat} containing protein accessions (default is \code{"Protein"}).}

\item{uniprot.output}{Character string specifying the type of output desired from UniProt, such as \code{"Genes"} (default) or any other valid UniProt field.}

\item{uniprot.extra}{Optional additional UniProt fields to retrieve (default is \code{NULL}).}

\item{uniprot.KEY}{The key type used for UniProt queries (default is \code{"UNIPROTKB"}).}

\item{uniprot.species.name}{Optional full species name used for UniProt queries (e.g., \code{"Homo sapiens"}). If not provided, \code{uniprot.taxId} must be given.}

\item{uniprot.taxId}{Optional UniProt taxonomy ID for the species of interest.}

\item{n.call}{Number of UniProt queries to make per batch (default is 99).}

\item{uniprot.file}{Path to a saved local annotation file (default is \code{"AnnotationsTable.xlsx"}).}

\item{shorten.gene.names}{Logical indicating whether long gene names should be shortened (default is \code{TRUE}).}

\item{split.char}{Character used to split multiple accessions (default is \code{";"}).}

\item{save.annotation}{Logical indicating whether to save the retrieved annotations to a file (default is \code{TRUE}).}

\item{save.annotation.path}{Path to the directory where the annotation file should be saved.}

\item{Data.path}{Path to the directory where data files are stored.}
}
\value{
A data table containing the original data with an additional \code{Genes} column with the annotated gene names.
}
\description{
This function annotates protein data by fetching gene information either from a local file or by querying the UniProt API.
The function is designed to handle protein accessions, split multi-accession inputs, and resolve gene names, optionally
shortening long gene names for easier visualization.
}
\details{
The function handles different annotation sources:
\itemize{
\item \code{"uniprot.file"}: Loads annotations from a previously saved Excel file.
\item \code{"annotation.file"}: Reads annotations from a user-provided local file.
\item \code{"uniprot"}: Queries UniProt directly to fetch annotations.
}

If the annotation source is \code{"uniprot"}, either \code{uniprot.taxId} or \code{uniprot.species.name} must be provided to ensure accurate species-specific annotation.

The function can handle inputs where multiple accessions are present in a single cell (proteins groups), typically separated by \code{";"}, and can split them for more precise annotations.

Gene names can be optionally shortened for better readability, especially for visualizations.
}
