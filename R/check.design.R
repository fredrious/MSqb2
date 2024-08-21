#' Check Design of Experimental Data
#'
#' This function verifies the design of experimental data provided in a data table. It checks for issues such as duplicate column names, multiple condition levels assigned to the same sample, and consistency in block and technical replicate assignments. Based on these checks, it determines the type of experimental design (e.g., fixed block, random block, or group comparison).
#'
#' @param dat A data.table containing experimental design metadata. It should include columns for `SampleID`, `Condition`, and optionally `Block` and `TechRep`.
#' @return A character string describing the type of experimental design detected and how the relevant factors are going to be handeled in the statistical model. Possible values include:
#' \itemize{
#'   \item "fixBlock.fixTechRep" - block: fixed-effect, TechRep: fixed-effect
#'   \item "fixBlock.randomTechRep" - block: fixed-effect, TechRep: random-effect
#'   \item "randomBlock.avgTechRep" - block: random-effect, average over TechRep.
#'   \item "fixBlock" - block: fixed-effect.
#'   \item "randomBlock" - block: random-effect.
#'   \item "fixTech" - TechRep: fixed-effect, no blocking.
#'   \item "randomTech" - TechRep: random-effect, no blocking.
#'   \item "group.comparison" - Simple group comparison design.
#' }
#' @examples
#' # Example data.table
#' dat <- data.table(SampleID = c("S1", "S2", "S3", "S4"),
#'                   Condition = c("C1", "C1", "C2", "C2"),
#'                   Block = c("B1", "B1", "B2", "B2"),
#'                   TechRep = c("T1", "T2", "T1", "T2"))
#' check.design(dat)
#'
#' @export
check.design <- function(dat) {

  if (any(duplicated(names(dat)))) {
    MSqb2:::.logg(
      level = FATAL,
      glue("There are duplicate column names in the metadata!")
    )
  }


  nCond.perSmp <- dat[, .("condition.levels" = uniqueN(Condition)), by = SampleID]
  if (any(nCond.perSmp$condition.levels > 1)) {

    # which samples match with more than 1 condition level
    dat[SampleID %in% nCond.perSmp[condition.levels > 1, SampleID]]
    MSqb2:::.logg(
      level = FATAL,
      glue("Multiple condition-levels assigned to the same SampleID!")
    )
  }


  if ("Block" %in% names(dat) & "TechRep" %in% names(dat)) {
    design.type <- .rand.fix.par(dat)
  } else {
    if ("Block" %in% names(dat)) {
      nCond.perSmpBlk <- dat[, .(nCondLvl = length(Condition)), by = c("SampleID", "Block")]
      Conds.perSmpBlk <- dat[, .(Conds.perSmp = toString(Condition)), by = c("SampleID", "Block")]
      is1to1 <- ifelse(any(nCond.perSmpBlk$nCondLvl > 1), FALSE, TRUE) # if 1to1 mapping between conditions and samples!

      if (is1to1) {
        design.type <- .check.paired.func(dat)
      } else {
        smp.lvl <- dcast(dat, SampleID ~ Condition, value.var = "Condition", fun.aggregate = length) %>%
          .[, !1] %>%
          unlist() %>%
          unique() %>%
          .[. > 0]
        dat[, TechRep := paste0("TechRep.", seq_along(Condition)), by = "SampleID"]
        dat[SampleID %in% {
          names(which(table(dat$SampleID) > 1))
        },
        SampleID := paste0(SampleID, ".", seq_along(Condition)),
        by = "SampleID"
        ]

        design.type <- .rand.fix.par(dat)
      }
    } else {
      if (!"TechRep" %in% names(dat)) {
        nCond.perSmp <- dat[, .(nCondLvl = length(Condition)), by = c("SampleID")]
        Conds.perSmp <- dat[, .(Conds.perSmp = toString(Condition)), by = c("SampleID")]
        is1to1 <- ifelse(any(nCond.perSmp[, !1] > 1), FALSE, TRUE) # if 1to1 mapping between conditions and samples!

        if (is1to1) {
          design.type <- "group.comparison"
        } else {
          smp.lvl <- dcast(dat, SampleID ~ Condition, value.var = "Condition", fun.aggregate = length) %>%
            .[, !1] %>%
            unlist() %>%
            unique() %>%
            .[. > 0]
          dat[, Block := SampleID]
          dat[, SampleID := paste0("Sample.", .I)]

          if (length(smp.lvl) == 1) {
            BlockSize <- uniqueN(dat$SampleID) / uniqueN(dat$Block)
            blk.cnd.smp <- dcast(dat, Block ~ Condition, value.var = "SampleID", \(x) paste(x, collapse = ","))
            design.type <- "fixBlock"
          } else {
            design.type <- "randomBlock"
            dat[, TechRep := paste0("TechRep.", seq_along(SampleID)), by = c("Block", "Condition")]
            if (uniqueN(dat$TechRep) == 1) dat[, TechRep := NULL]
          }
        }
      } else {
        nCond.perSmpTech <- dat[, .(nCondLvl = length(Condition)), by = c("SampleID", "TechRep")]
        Conds.perSmpTech <- dat[, .(Conds.perSmp = toString(Condition)), by = c("SampleID", "TechRep")]

        if (any(nCond.perSmpTech$nCondLvl > 1)) { # blocking in data

          dat[, Block := SampleID]
          dat[, SampleID := paste0("Sample.", .I)]
          design.type <- .rand.fix.par(dat)
        } else { # no blocking

          design.type <- .rand.fix.par(dat)
        }
      }
    }
  }
  return(design.type)
}




.rand.fix.par <- function(dat) {
  if ("TechRep" %in% names(dat)) {
    Conds.perTech <- dat[, .(nCondLvl = length(Condition)), by = c("TechRep")]
    isSymmetricTech <- ifelse(uniqueN(Conds.perTech$nCondLvl) == 1, TRUE, FALSE)
  } else {
    isSymmetricTech <- NULL
  }

  if ("Block" %in% names(dat)) {
    Conds.perBlock <- dat[, .(nCondLvl = length(Condition)), by = c("Block")]
    isSymmetricBlock <- ifelse(uniqueN(Conds.perBlock$nCondLvl) == 1, TRUE, FALSE)
  } else {
    isSymmetricBlock <- NULL
  }


  if (!is.null(isSymmetricBlock) & !is.null(isSymmetricTech)) {
    if (isSymmetricBlock & isSymmetricTech) design.type <- "fixBlock.fixTechRep"
    if (isSymmetricBlock & !isSymmetricTech) design.type <- "fixBlock.randomTechRep"
    if (!isSymmetricBlock & !isSymmetricTech) design.type <- "randomBlock.avgTechRep"
  }

  if (!is.null(isSymmetricBlock) & is.null(isSymmetricTech)) {
    design.type <- ifelse(isSymmetricBlock, "fixBlock", "randomBlock")
  }

  if (is.null(isSymmetricBlock) & !is.null(isSymmetricTech)) {
    design.type <- ifelse(isSymmetricTech, "fixTech", "randomTech")
  }

  if (is.null(isSymmetricBlock) & is.null(isSymmetricTech)) {
    design.type <- "group.comparison"
  }

  return(design.type)
}



.check.paired.func <- function(dat) {
  blk.cnd.smp <- dcast(dat, Block ~ Condition, value.var = "SampleID", \(x) paste(x, collapse = ","))
  blk.cnd.cnd <- dcast(dat, Block ~ Condition, value.var = "Condition", \(x) paste(unique(x), collapse = ","))
  blk.cnd.cnt <- dcast(dat, Block ~ Condition, value.var = "Condition", fun.aggregate = length)

  split(blk.cnd.cnt, by = "Block", keep.by = FALSE) %>%
    Reduce(\(x, y) {
      if (identical(x, y)) x else FALSE
    }, .) %>%
    is.logical(.) %>%
    ifelse(., "randomBlock", "fixBlock") %>%
    return()
}
