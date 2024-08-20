

exportStatTable <- function(dat,
                           lfc.col = "logFC",
                           apv.col = "adj.P.Val",
                           gene.col = "Gene",
                           prot.col = "Protein",
                           negFill = "skyblue1",
                           posFill = "pink",
                           noSigFill = "grey90",
                           sigFill = "red",
                           sigFont = "darkred",
                           tableStyle = "TableStyleLight8",
                           file.name = "limma.stat_output",
                           path = Tables.path,
                           sig.level = 0.05,
                           silent = FALSE) {
  negStyle <- createStyle(fontColour = "black", bgFill = negFill)
  posStyle <- createStyle(fontColour = "black", bgFill = posFill)
  noSigstyle <- createStyle(fontColour = "black", bgFill = noSigFill)
  sigStyle.sub <- createStyle(fontColour = sigFont, textDecoration = "bold")
  sigStyle.rest <- createStyle(fontColour = "black", textDecoration = "bold")
  sigLevelStyle <- createStyle(border = "Top", borderColour = sigFont, borderStyle = "medium")
  headStyle <- createStyle(
    fontColour = "grey10", bgFill = "white", fgFill = "white",
    halign = "center", valign = "center", textDecoration = "Bold",
    border = "TopBottomLeftRight"
  )


  dat <- dat[order(Comparison, rank.id)]

  if (!silent) {
    message(
      "\n\nBy defualt, the statistical parameters corresponding to each comparison ",
      "will be saved in a separate sheet in the excel file:\n ",
      file.path(path, paste0(file.name, ".xlsx")),
      "\nComparison names will be used as sheet names. The maximum allowed length ",
      "for the sheet names is 30 charachters!"
    )
  }

  compdt <- dat[, .(Comparison = unique(as.character(Comparison)))]
  compdt[, nch := nchar(Comparison)]
  compdt[, Comp := ifelse(nch > 30, paste0("Comp", seq_along(Comparison)), Comparison)]


  if (any(compdt$nch > 30)) {
    if (!silent) {
      message(
        "\n\n!!! The following comparisons have names longer than 30 charachters! ",
        "To save the statistical parameters in the excel file, these comparisons will be renamed ",
        "according to the table below (also saved as txt file: 'longCompNames.txt'\n",
        "The original names will still be used for data visualization."
      )
    }

    print(compdt[nch > 30, .(Comparison, Comp)])
    fl <- file.path(path, "longCompNames.txt")
    sink(fl)
    cat("Comparisons with names longer than 30 charachters:\n==========\n")
    print(compdt[nch > 30, .(Comparison, Comp)])
    sink()
  }



  wb <- createWorkbook()

  for (ic in 1:nrow(compdt)) {
    comp <- compdt[ic, Comparison]
    sheet <- make.names(as.character(compdt[ic, Comp]))

    subdat <- dat[Comparison == comp]
    nrw <- nrow(subdat)

    lfc <- which(names(subdat) == lfc.col)
    apv <- which(names(subdat) == apv.col)
    prt <- which(names(subdat) == prot.col)
    gn <- which(names(subdat) == gene.col)

    min.apv <- min(subdat[, ..apv.col], na.rm = TRUE)
    sig.row.sub <- tail(subdat[get(apv.col) < sig.level, .I[sig.level == max(sig.level)]], 1) + 1
    sig.prot.sub <- as.character(unique(subdat[get(apv.col) < sig.level, ..prot.col]))
    sig.prot.rest <- setdiff(as.character(unique(dat[get(apv.col) < sig.level, ..prot.col])), sig.prot.sub)
    sig.row.rest <- which(unlist(subdat[, ..prot.col]) %in% sig.prot.rest)



    addWorksheet(wb, sheet)
    freezePane(wb, sheet = sheet, firstRow = TRUE) ## freeze first row and column

    writeDataTable(wb, sheet, subdat,
      colNames = TRUE,
      withFilter = FALSE,
      tableStyle = tableStyle,
      headerStyle = headStyle
    )


    # adj.p color formatting
    if (min.apv < sig.level) {
      conditionalFormatting(wb, sheet,
        cols = apv,
        rows = 2:nrw,
        style = c(sigFill, noSigFill),
        rule = c(min.apv, sig.level),
        type = "colourScale"
      )
    } else {
      conditionalFormatting(wb, sheet, cols = apv, rows = 2:nrw, rule = "<=1", style = noSigstyle)
    }


    # log fold change color formatting
    conditionalFormatting(wb, sheet, cols = lfc, rows = 2:nrw, rule = "<0", style = negStyle)
    conditionalFormatting(wb, sheet, cols = lfc, rows = 2:nrw, rule = ">0", style = posStyle)


    # highlight sig proteins in other comparisons! "black"
    if (length(sig.row.rest) > 0) {
      addStyle(wb, sheet, cols = prt, rows = sig.row.rest, style = sigStyle.rest)
      addStyle(wb, sheet, cols = gn, rows = sig.row.rest, style = sigStyle.rest)
    }


    # highlight sig proteins "darkred"
    if (!is.na(sig.row.sub)) {
      addStyle(wb, sheet, cols = prt, rows = 2:sig.row.sub, style = sigStyle.sub)
      addStyle(wb, sheet, cols = gn, rows = 2:sig.row.sub, style = sigStyle.sub)
      ## horizontal line marking the sig. level
      addStyle(wb, sheet, sigLevelStyle, rows = sig.row.sub + 1, cols = 1:ncol(dat), gridExpand = TRUE)
    }
  }

  saveWorkbook(wb, file = file.path(path, paste0(file.name, ".xlsx")), overwrite = TRUE)
  paste0("Statical parameters were saved in file: ", file.path(path, paste0(file.name, ".xlsx")))

  return(wb)
}
