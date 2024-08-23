#' Merge Metadata with Mass Spectrometry Data
#'
#' This function merges metadata with mass spectrometry data, ensuring consistency between the pools and channels in both datasets.
#' It performs checks to validate the presence and consistency of `Pool` and `Channel` information between the datasets, and adjusts
#' them as necessary to match.
#'
#' @param msDT A `data.table` containing the mass spectrometry data. This table should include columns such as `Pool`, `Channel`, and `Protein`.
#' @param metadata A `data.table` containing the metadata to be merged with `msDT`. This table should include columns such as `Pool` and `Channel`.
#' @param extMixFrac A `data.table` containing additional information for external mix fractions, used if `Pool` information is missing in the metadata.
#' @param int.nm A character string specifying the internal name used for processing. (This parameter is referenced in the function but not used within the provided code.)
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item `"msDT"`: The merged mass spectrometry data as a `data.table`.
#'   \item `"metadata"`: The merged metadata as a `data.table`.
#' }
#'
#' @details The `merge_ms_metadata` function ensures that the `Pool` and `Channel` columns in the mass spectrometry data (`msDT`) and the metadata are consistent.
#' It checks for the presence of `Pool` information in both datasets, adjusts levels if necessary, and merges the metadata into the mass spectrometry data.
#'
#' Specific steps include:
#' \itemize{
#'   \item Verifying and correcting `Pool` levels between `msDT` and `metadata`.
#'   \item Handling cases where `Pool` information is missing in either `msDT` or `metadata`.
#'   \item Ensuring consistency of `Channel` information between the datasets.
#'   \item Merging the metadata into `msDT` with consideration for potential inconsistencies.
#' }
#'
#' @importFrom magrittr %>%
#' @export
merge_ms_metadata <- function(msDT, metadata, extMixFrac, int.nm) {
  char2fact(metadata)
  char2fact(msDT)

  ## CHECK THE POOL
  ## ---------------

  # check if the pools are the same in PSM data and in dataset
  if ("Pool" %in% names(msDT)) {
    poollvl <- msDT$Pool %>%
      unique() %>%
      unlist()
  }

  # check if the pools are the same in PSM data and in metadata
  if ("Pool" %in% names(metadata)) {
    poollvl.exp <- metadata$Pool %>%
      unique() %>%
      unlist()
  }




  # if pool not exists in design, but in data
  if ("Pool" %in% names(metadata) & !"Pool" %in% names(msDT)) {
    if (length(poollvl.exp) > 1) {
      print(metadata)
      stop(
        "\nThe Pool column in the metadata has multiple levels while dataset does not contain Pool information!",
        "\nPool levels in the metadata: ", paste(as.character(poollvl.exp), collapse = ", "),
        "\nPlease edit the metadata file and re-run the workflow."
      )
    }
  }


  # if pool not exists in design, but in data
  if (!"Pool" %in% names(metadata) & "Pool" %in% names(msDT)) {
    if (length(poollvl) > 1) {
      print(metadata)
      stop(
        "\nThe Pool column in the dataset has multiple levels while the metadata does not contain Pool information!",
        "\nPool levels in the dataset: ", paste(as.character(poollvl), collapse = ", "),
        "\nPlease edit the metadata file and re-run the workflow."
      )
    } else { # add pool to metadata from extMixFrac
      warning("\nThe Pool column does not exist in metadata. It will be added from the input dataset.")
      metadata[, Pool := unique(extMixFrac$Pool)]
      poollvl.exp <- unique(extMixFrac$Pool)
    }
  }


  # if pool exists in both
  if ("Pool" %in% names(metadata) & "Pool" %in% names(msDT)) {
    if (all(!poollvl %in% poollvl.exp)) {
      for (ip in poollvl[!poollvl %in% poollvl.exp]) { # correcting the pool levels in exp design
        ipx <- menu(sort(poollvl.exp),
          title = paste(
            "\nPool levels in the dataset do not match the Pool levels in the metadata!",
            "\nPool levels in the dataset: ", paste(as.character(poollvl), collapse = ", "),
            "\nPool levels in the metadata: ", paste(as.character(poollvl.exp), collapse = ", "),
            "\nWhich Pool in the metadata corresponds to", ip, "in the dataset?"
          )
        )
        msDT[Pool == ip, Pool := poollvl.exp[ipx]]
        poollvl <- msDT$Pool %>%
          unique() %>%
          unlist()
      }
    }

    if (length(poollvl) > length(poollvl.exp)) { # if the pool levels in the data are more than those in Exp design
      ipx <- menu(c("yes", "no"),
        title = paste(
          "\nThe Pool levels in the dataset are more than the Pool levels in the metadata!",
          "\nPool levels in the dataset: ", paste(as.character(poollvl), collapse = ", "),
          "\nPool levels in the metadata: ", paste(as.character(poollvl.exp), collapse = ", "),
          "\nSubset the dataset by excluding the rows corresponding to pool levels",
          setdiff(poollvl.exp, poollvl), "? \n(by chosing no the workflow will break!)"
        )
      )
      if (ipx == 1) {
        msDT <- msDT[Pool %in% poollvl.exp] # subset data
      } else {
        stop()
      }
    }


    if (length(poollvl) < length(poollvl.exp)) { # if the pool levels in the Exp design are more than those in dataset
      message(
        "\nThe Pool levels in the metadata are more than the Pool levels in the dataset!",
        "\nPool levels in the dataset: ", paste(as.character(poollvl), collapse = ", "),
        "\nPool levels in the metadata: ", paste(as.character(poollvl.exp), collapse = ", "),
        "\nBy default, the extra levels in the metadata will be ignored!"
      )
    }
  }




  ## CHECK THE CHANNEL
  ## ------------------
  metadata[, Channel := as.factor(as.character(Channel))]
  msDT[, Channel := as.factor(as.character(Channel))]

  x <- as.character(sort(unique(metadata$Channel)))[1]
  y <- as.character(sort(unique(msDT$Channel)))
  which(y %like% x) %>% y[.] -> y
  l <- c(x, y)
  exchar <- strsplit2(
    names(which.max(sapply(l, nchar))),
    names(which.min(sapply(l, nchar)))
  )

  if (nchar(exchar) > 0) {
    if (y %like% exchar) {
      msDT[, Channel := gsub(exchar, "", Channel, fixed = TRUE)]
    } else {
      metadata[, Channel := gsub(exchar, "", Channel, fixed = TRUE)]
    }
  }

  metadata[, Channel := as.factor(as.character(Channel))]
  char2fact(metadata)
  msDT[, Channel := as.factor(as.character(Channel))]
  char2fact(msDT)


  if (any(!names(metadata) %in% names(msDT))) {
    msDT <- merge(msDT, metadata, all.y = TRUE)
  }

  # above: all.y = TRUE set because some samples might need to be removed from data.
  # But it may also lead to extra NA in other columns. These must be excluded.
  msDT <- msDT[!is.na(Protein), ]

  return(list("msDT" = msDT, "metadata" = metadata))
}
