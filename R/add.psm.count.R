#' Add PSM count and related information to msDT
#'
#' This function calculates and adds PSM counts to the msDT 
#'
#' @param dt A msDT containing the original data.
#' @param param The parameter or column name used for grouping.
#'
#' @return A modified data.table with added columns representing PSM count-related information.
#'
#' @examples
add.psm.count <- function(dt, param) {
  # psm.df <- dt[, .(PSM.count = uniqueN(Feature)), by = param][order(PSM.count)]
  psm.df <- dt[, c(..param, "PSM.count")][order(PSM.count)] %>% unique()
  
  # psm.df <-  psm.df[order(PSM.count), .(PSM.count)] %>% unique() %>% 
  #   .[, PSM.interval := shift(PSM.count, fill = PSM.count[.N]+1, type = "lead") - PSM.count] %>% 
  #   merge(psm.df, ., by = "PSM.count")
  
  
  psm.df <- dt[, .(Abundance.mean = mean(Abundance, na.rm = T), PSM.count = PSM.count), by = param] %>% 
    unique() %>% 
    .[order(PSM.count, Abundance.mean)] %>% 
    .[, PSM.mean := (seq_along(get(param))/.N)+PSM.count, by = PSM.count]
  
  dt <- merge(dt, psm.df[, c(param, "PSM.mean"), with = FALSE], by = "Protein")
  return(dt)
}
