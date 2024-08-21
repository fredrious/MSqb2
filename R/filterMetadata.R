#' Filters metadata.
#'
#' If specific levels or combination of levels are to be filtered from the input data data.
#' In 'filters', if a list-element is a vector of length one, the level will be completely removed
#' from the input data. If list-elements have length larger than one, the combination of the
#' levels in the vector will be used to filter data. Column names do not need to be passed,
#' as column-matching will be done automatically. See examples for more details.
#'
#' @param dat data.table or data.frame of the categorical variables (metadata).
#' @param filters list of character vectors containing levels in the input data.
#'
#' @return filtered data
#' @export
#'
#' @examples
#' dat <- data.frame(Condition = rep(rep(paste0("Cond.", 1:2), each = 5), 2) ,
#'                  BioRep = rep(rep(paste0("Rep.", 1:5)), 2),
#'                  Pool = rep(paste0("Pool.", 1:2), each = 10) )
#'
#' # filter Rep.2 OR Pool.1 completely (both leveles will be removed).
#' filter <- list("Rep.2", "Pool.1")
#' fdt <- filterMetadata(dat, filter)
#'
#'
#' # filter rows matched by: Rep.5 AND Cond.2
#' filter <- list(c("Rep.5", "Cond.2"))
#' fdt <- filterMetadata(dat, filter)
#'
#'
#' # filter rows matched by: Pool.1 AND Cond.1 AND Rep.3
#' filter <- list(c("Pool.1", "Cond.1", "Rep.3"))
#' fdt <- filterMetadata(dat, filter)
#'
#'
#' # filter rows matched by: Pool.1 AND Cond.1 AND Rep.3 OR Pool.1 AND Cond.1 AND Rep.4
#' filter <- list(c("Pool.1", "Cond.1", "Rep.3"), c("Pool.1", "Cond.1", "Rep.4"))
#' fdt <- filterMetadata(dat, filter)
#'
#'
#' # filter rows matched by: Cond.1 AND Rep.3 OR Cond.1 AND Rep.5
#' filter <- list(c("Cond.1", "Rep.3"), c("Cond.1", "Rep.5"))
#' fdt <- filterMetadata(dat, filter)
#'
#'
#' # filter rows matched by: Rep.5 AND Cond.2 OR Pool.1 AND Cond.1
#' filter <- list(c("Rep.5", "Cond.2"), c("Pool.1", "Cond.1"))
#' fdt <- filterMetadata(dat, filter)


filterMetadata <- function(dat, filters) {
  dat %<>% as.data.table(.)
  dat <- fact2char(dat)
  if (!is.null(filters)) {
    cond <- NULL
    for (p in seq_along(filters)) {
      cnd <- NULL
      for (ip in filters[[p]]) {
        cl <- names(dat)[dat[, .SD %like% ip]]
        if (length(cl) > 0) {
          cnd <- c(cnd, paste0(cl, " == '", ip, "'"))
        }
      }
      cnd <- paste(cnd, collapse = " & ")
      cond <- c(cond, cnd)
    }
    if (length(cond) > 0) {
      cond <- paste(cond, collapse = " | ")
    } else {
      cond <- cnd
    }
    dat[!(eval(parse(text = cond))), ] %>%
      MSqb2::char2fact(.) %>%
      return(.)
  } else {
    MSqb2::char2fact(dat) %>%
      return(.)
  }
}

