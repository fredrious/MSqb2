
################################
## function to i. calculate PSM ratios
##             ii. impute NA for PSM within Peptide, Pool, Run and Channel
##             iii. return original Abundance, imputed Abundance & median Abindance
## note: median abundance inserted in datatable for both charge states, meaning median abundance is duplicated for each peptide.
## This is useful for ggplot to have all three levels of Abundance in the same datatable.

f.psm2ppt <- function(ww) {
  ww <- t(ww[, "Feature"])
  colnames(ww) <- t(ww[, "Feature"])
  # filling NAs of 2nd Charge by nonNAs of 1st Charge
  ww[is.na(ww[, 2]), 2] <- ww[is.na(ww[, 2]), 1] - median(ww[, 1] - ww[, 2], na.rm = TRUE)
  # filling NAs of 1st Charge by nonNAs of 2nd Charge
  ww[is.na(ww[, 1]), 1] <- ww[is.na(ww[, 1]), 2] + median(ww[, 1] - ww[, 2], na.rm = TRUE)
  ww.med <- data.table(nm = rownames(ww), ww, medPeptide = apply(ww, 1, function(x) median(x, na.rm = TRUE))) # table of imputed and mean PSM

  PSM.imp <- melt.data.table(ww.med[, -4],
    id.vars = "nm",
    variable.name = "Feature", value.name = "impAbundance"
  )
  PSM.imp[, c("Pool", "Run", "Channel") := tstrsplit(nm, "_", fixed = TRUE)] # table of imputed abundance
  PSM.imp[, c("Protein", "Peptide", "Charge") := tstrsplit(Feature, "_", fixed = TRUE)]

  PSM.med <- melt.data.table(ww.med[, c(1, 4)], id.vars = "nm", value.name = "medAbundance")
  PPT.dt <- PSM.imp[PSM.med[, -"variable"], on = "nm"] # add median abundance

  return(PPT.dt)
}
################################
