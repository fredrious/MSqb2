StatTable <- function(fit_eb,
                      dt,
                      featureData = NULL, # featuredata, if provided
                      variable = "Protein",
                      Gene.name = "Genes",
                      adjust = "BH",
                      
                      # split.char = ";",
                      save2file = TRUE,
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
    tab.full <- topTable(fit_eb, coef = j, number = Inf, adjust = adjust, confint = TRUE, genelist = NULL)
    keepnm <- ifelse(variable %in% names(tab.full), FALSE, variable)
    tab.full <- data.table(tab.full, keep.rownames = keepnm)
    
    # to add comparison and id column at the beginning of the datatable
    tab.full <- data.table(Comparison = colnames(fit_eb$contrasts)[j], tab.full)
    tab.full <- tab.full[order(Comparison, P.Value, logFC, decreasing = FALSE), ]
    tab.full <- data.table(rank.id = seq(nrow(tab.full)), tab.full)
    topList <- rbind(topList, tab.full)
  }
  topList <- topList[!is.na(Comparison), ]
  MSqb2::char2fact(topList)
  
  ## merge toplist with dt
  if (!is.null(dt) && length(dt) != 0) {
    if (!Gene.name %in% names(topList)) {
      topList <- names(dt) %>%
        match(c(variable, Gene.name, "Peptide.count", "PSM.count"), .) %>%
        .[!is.na(.)] %>%
        unique(.) %>%
        dt[, .SD, .SDcols = .] %>%
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
  
  
  
  if (save2file) {
    try(exportStatTable(
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
