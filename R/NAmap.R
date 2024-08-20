

NAmap <- function(dt,
                  val = "Abundance",
                  row.p = "Protein",
                  col.p = "auto",
                  export.as = "long", # "long", "wide", "matrix"
                  sp = "_.._") {


  # subset available categorical parameters
  if (is.null(col.p) || tolower(col.p) == "auto") {
    col.p <- intersect(
      unique(c("Filename", "SampleID", "Pool", "Condition", "Channel", "BioRep", "TechRep", "Fraction", col.p)),
      names(dt))
  }

  
  # #order data
  setorderv(dt, cols = col.p)

  # generate heatmap annotation
  ha <- unique(dt[, ..col.p])
  ha[, cls := do.call(paste, c(.SD, sep = sp)), .SDcols = col.p]

  # convert to wide and then to Matrix
  dna <- dt[, .(Absm = sum(get(val))), by = c(row.p, col.p)]
  form <- paste0(
    ifelse(length(row.p) > 1, paste(row.p, collapse = "+"), row.p),
    "~",
    ifelse(length(col.p) > 1, paste(col.p, collapse = "+"), col.p)
  )

  dna <- dcast.data.table(dna,
    formula = form,
    value.var = "Absm",
    sep = sp
  ) # cast input data -- wide format

  setcolorder(dna, c(row.p, ha$cls))
  rowtg <- dna[, ..row.p]
  dna <- dna[, lapply(.SD, function(x) !is.finite(x)), .SDcol = setdiff(names(dna), row.p)]
  dna <- cbind(rowtg, dna)

  ixNA <- which(!apply(dna[, -..row.p], 1, function(x) any(!x)))
  if (length(ixNA) > 0) {
    message(paste0("Following ", row.p, "s have complete missingness and will be removed.\n"))
    xNA <- dna[ixNA, ..row.p]
    print(xNA)
    dna <- dna[setdiff(seq(nrow(dna)), ixNA), ]
  }

  if (toupper(export.as) == "MATRIX") dna <- as.matrix(data.frame(dna, row.names = row.p))

  if (toupper(export.as) == "LONG") {
    dna <- MSqb2:::.w2l(dna, col.p = col.p, val = val, row.p = row.p, sp = sp)
    setnames(dna, val, "isNA")
  }

  return(dna)
}
