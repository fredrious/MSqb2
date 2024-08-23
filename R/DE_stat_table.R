#' Generate result table of the Differential Expression Analysis
#'
#' This function generates a differential expression statistics table from a fitted empirical Bayes model, merges it with additional data if provided, and optionally saves the results to a file.
#'
#' @param fit_eb A fitted empirical Bayes model object (from `limma`).
#' @param dat A data.table containing additional data to merge with the statistics table.
#' @param featureData Optional data.table with feature data (gene names, ...) to be merged with the statistics table.
#' @param variable A character string specifying the variable name for proteins or features in the data. Default is "Protein".
#' @param Gene.name A character string specifying the column name for gene names. Default is "Genes".
#' @param adjust A character string specifying the method for p-value adjustment. Default is "BH".
#' @param save Logical, if `TRUE`, saves the result to a file. Default is `TRUE`.
#' @param path A character string specifying the path to save the file. Default is `Tables.path`.
#' @param file.name A character string specifying the file name for saving the output. Default is "StatTest_output".
#' @param ... Additional arguments passed to the `export_DE_result` function.
#'
#' @return A data.table containing the differential expression statistics.
#'
#' @importFrom magrittr %>% %<>%
#' @export

DE_stat_table <- function(fit_eb,
                          dat,
                          featureData = NULL, # featuredata, if provided
                          variable = "Protein",
                          Gene.name = "Genes",
                          adjust = "BH",

                          # split.char = ";",
                          save = TRUE,
                          path = Tables.path,
                          file.name = "StatTest_output",
                          ...) {

  colsort <- c(
    "rank.id", "Comparison", Gene.name, variable, "Peptide", "logFC",
    "P.Value", "adj.P.Val", "AveExpr", "Conditions"
  )



  topList <- NULL
  for (j in 1:ncol(fit_eb$contrasts)) {
    ## topTable
    tab.full <- topTable(fit_eb, coef = j, number = Inf, adjust.method = adjust, confint = TRUE, genelist = NULL)
    keepnm <- ifelse(variable %in% names(tab.full), FALSE, variable)
    tab.full <- data.table(tab.full, keep.rownames = keepnm)

    # to add comparison and id column at the beginning of the datatable
    tab.full <- data.table(Comparison = colnames(fit_eb$contrasts)[j], tab.full)
    tab.full <- tab.full[order(Comparison, P.Value, logFC, decreasing = FALSE), ]
    tab.full <- data.table(rank.id = seq(nrow(tab.full)), tab.full)
    topList <- rbind(topList, tab.full)
  }
  topList <- topList[!is.na(Comparison), ]
  char2fact(topList)

  ## merge toplist with dat
  if (!is.null(dat) && length(dat) != 0) {
    if (!Gene.name %in% names(topList)) {
      topList <- names(dat) %>%
        match(c(variable, Gene.name, "Peptide.count", "PSM.count"), .) %>%
        .[!is.na(.)] %>%
        unique(.) %>%
        dat[, .SD, .SDcols = .] %>%
        unique(.) %>%
        merge.data.table(topList, ., by = variable)
    }
  } else {
    colsort <- setdiff(colsort, Gene.name)
    Gene.name <- NULL
  }



  if ("Peptide.count" %in% names(topList)) colsort <- c(colsort, "Peptide.count", "PSM.count")

  topList.s <- names(topList) %>%
    match(., colsort) %>%
    .[!is.na(.)] %>%
    sort(.) %>%
    colsort[.] %>%
    topList[, .SD, .SDcols = .]


  # if featureData provided, merge
  if (!is.null(featureData)) {
    topList.s %<>% merge.data.table(., featureData,
                                    by = variable,
                                    all.x = TRUE,
                                    all.y = FALSE,
                                    sort = FALSE)
  }


  # include condition names in the statlist
  # merge with condition levels
  ContComp <- fit_eb$contrasts %>%
    t(.) %>%
    apply(., 1, \(x) paste( unique(names(x[x!=0])), collapse = ";") ) %>%
    data.frame(Conditions = .) %>%
    as.data.table(., keep.rownames = "Comparison")
  topList.s %<>% merge(., ContComp, by = "Comparison")



  if (save) {
    try(export_DE_result(
      dat = topList.s[order(Comparison, rank.id), ],
      path = path,
      file.name = file.name,
      gene.col = Gene.name,
      prot.col = variable, ...
    ))
  }
  rm(topList.s)

  # topList <- topList[, -"Comparison"]
  # setnames(topList, old = "Comparison.org", new = "Comparison")

  topList %>% merge.data.table(., ContComp, by = "Comparison") %>%
    .[order(Comparison, rank.id)] %>%
    return()
}
