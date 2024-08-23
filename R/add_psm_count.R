#' Add PSM Count and Related Information to msDT
#'
#' This function calculates and adds PSM (Peptide-Spectrum Match) counts to the provided mass spectrometry data table (`msDT`).
#'
#' @param dat A `data.table` containing the original MS data, typically with columns like `Feature`, `Protein`, `Abundance`, etc.
#' @param param A character string specifying the column name used for grouping when calculating PSM counts (e.g., `"Protein"`).
#'
#' @return A modified `data.table` with additional columns representing PSM count-related information, including mean PSM counts and abundance statistics.
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item Groups the data by the specified `param` and calculates the PSM count.
#'   \item Computes the mean abundance for each group.
#'   \item Adds a column representing the mean PSM count adjusted by the abundance rank.
#'   \item Merges this information back into the original data table (`msDT`).
#' }
#'
#' @importFrom magrittr %>%
#' @export
add_psm_count <- function(dat, param) {
  psm.df <- dat[, c(..param, "PSM.count")][order(PSM.count)] %>% unique()

  psm.df <- dat[, .(Abundance.mean = mean(Abundance, na.rm = T), PSM.count = PSM.count), by = param] %>%
    unique() %>%
    .[order(PSM.count, Abundance.mean)] %>%
    .[, PSM.mean := (seq_along(get(param))/.N)+PSM.count, by = PSM.count]

  dat <- merge(dat, psm.df[, c(param, "PSM.mean"), with = FALSE], by = "Protein")
  return(dat)
}
