

check.design <- function(pdt) {
  if (any(duplicated(names(pdt)))) {
    MSqb2:::.logg(
      level = FATAL,
      glue("There are duplicate column names in the metadata!")
    )
  }


  nCond.perSmp <- pdt[, .("condition.levels" = uniqueN(Condition)), by = SampleID]
  if (any(nCond.perSmp$condition.levels > 1)) {

    # which samples match with more than 1 condition level
    pdt[SampleID %in% nCond.perSmp[condition.levels > 1, SampleID]]
    MSqb2:::.logg(
      level = FATAL,
      glue("Multiple condition-levels assigned to the same SampleID!")
    )
  }


  if ("Block" %in% names(pdt) & "TechRep" %in% names(pdt)) {
    design.type <- .rand.fix.par(pdt)
  } else {
    if ("Block" %in% names(pdt)) {
      nCond.perSmpBlk <- pdt[, .(nCondLvl = length(Condition)), by = c("SampleID", "Block")]
      Conds.perSmpBlk <- pdt[, .(Conds.perSmp = toString(Condition)), by = c("SampleID", "Block")]
      is1to1 <- ifelse(any(nCond.perSmpBlk[, !1] > 1), FALSE, TRUE) # if 1to1 mapping between conditions and samples!

      if (is1to1) {
        design.type <- .check.paired.func(pdt)
      } else {
        smp.lvl <- dcast(pdt, SampleID ~ Condition, value.var = "Condition", fun.aggregate = length) %>%
          .[, !1] %>%
          unlist() %>%
          unique() %>%
          .[. > 0]
        pdt[, TechRep := paste0("TechRep.", seq_along(Condition)), by = "SampleID"]
        pdt[SampleID %in% {
          names(which(table(pdt$SampleID) > 1))
        },
        SampleID := paste0(SampleID, ".", seq_along(Condition)),
        by = "SampleID"
        ]

        design.type <- .rand.fix.par(pdt)
      }
    } else {
      if (!"TechRep" %in% names(pdt)) {
        nCond.perSmp <- pdt[, .(nCondLvl = length(Condition)), by = c("SampleID")]
        Conds.perSmp <- pdt[, .(Conds.perSmp = toString(Condition)), by = c("SampleID")]
        is1to1 <- ifelse(any(nCond.perSmp[, !1] > 1), FALSE, TRUE) # if 1to1 mapping between conditions and samples!

        if (is1to1) {
          design.type <- "group.comparison"
        } else {
          smp.lvl <- dcast(pdt, SampleID ~ Condition, value.var = "Condition", fun.aggregate = length) %>%
            .[, !1] %>%
            unlist() %>%
            unique() %>%
            .[. > 0]
          pdt[, Block := SampleID]
          pdt[, SampleID := paste0("Sample.", .I)]

          if (length(smp.lvl) == 1) {
            BlockSize <- uniqueN(pdt$SampleID) / uniqueN(pdt$Block)
            blk.cnd.smp <- dcast(pdt, Block ~ Condition, value.var = "SampleID", \(x) paste(x, collapse = ","))
            design.type <- "fixBlock"
          } else {
            design.type <- "randomBlock"
            pdt[, TechRep := paste0("TechRep.", seq_along(SampleID)), by = c("Block", "Condition")]
            if (uniqueN(pdt$TechRep) == 1) pdt[, TechRep := NULL]
          }
        }
      } else {
        nCond.perSmpTech <- pdt[, .(nCondLvl = length(Condition)), by = c("SampleID", "TechRep")]
        Conds.perSmpTech <- pdt[, .(Conds.perSmp = toString(Condition)), by = c("SampleID", "TechRep")]

        if (any(nCond.perSmpTech$nCondLvl > 1)) { # blocking in data

          pdt[, Block := SampleID]
          pdt[, SampleID := paste0("Sample.", .I)]
          design.type <- .rand.fix.par(pdt)
        } else { # no blocking

          design.type <- .rand.fix.par(pdt)
        }
      }
    }
  }
  return(design.type)
}




.rand.fix.par <- function(pdt) {
  if ("TechRep" %in% names(pdt)) {
    Conds.perTech <- pdt[, .(nCondLvl = length(Condition)), by = c("TechRep")]
    isSymmetricTech <- ifelse(uniqueN(Conds.perTech$nCondLvl) == 1, TRUE, FALSE)
  } else {
    isSymmetricTech <- NULL
  }

  if ("Block" %in% names(pdt)) {
    Conds.perBlock <- pdt[, .(nCondLvl = length(Condition)), by = c("Block")]
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



.check.paired.func <- function(pdt) {
  blk.cnd.smp <- dcast(pdt, Block ~ Condition, value.var = "SampleID", \(x) paste(x, collapse = ","))
  blk.cnd.cnd <- dcast(pdt, Block ~ Condition, value.var = "Condition", \(x) paste(unique(x), collapse = ","))
  blk.cnd.cnt <- dcast(pdt, Block ~ Condition, value.var = "Condition", fun.aggregate = length)

  split(blk.cnd.cnt, by = "Block", keep.by = FALSE) %>%
    Reduce(\(x, y) {
      if (identical(x, y)) x else FALSE
    }, .) %>%
    is.logical(.) %>%
    ifelse(., "randomBlock", "fixBlock") %>%
    return()
}
