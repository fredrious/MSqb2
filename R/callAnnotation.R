#' call Annotation using AnnotationDbi
#'
#' @param prtOrganism AnnotationDb object (genome). Currently one of human, mouse and yeast
#' @param annot.keys keys to select records for from the database, such as ENSEMBL ids.
#' @param keytype keytype that matches the keys used. For instance, if annot.keys are ensml ids, the keytype must be ENSEMBL"
#' @param annot.cols the columns or kinds of things that can be retrieved from the database. "SYMBOL","GENENAME","ENTREZID", "REFSEQ", ...
#' @param AnnDbi.multiVals What should mapIds do when there are multiple values that could be returned? look at ?AnnotationDbi::mapIds > multiVals
#'
#' @return a data.table with annot.keys and annot.cols
#' @export
#'
callAnnotationDbi <- function(prtOrganism,
                               annot.keys = gene_names,
                               keytype = "ENSEMBL",
                               annot.cols = c("SYMBOL","GENENAME","ENTREZID", "REFSEQ"),
                               AnnDbi.multiVals = "list") {
   org <- toupper(prtOrganism)
   if (org == "YEAST") {
      annot.cols = c("COMMON", "GENENAME", "ENTREZID", "REFSEQ")
   }

   ## prepare genome
   genome_data <- switch(
      org,
      "HUMAN" = {require(org.Hs.eg.db); org.Hs.eg.db},
      "MOUSE" = {require(org.Mm.eg.db); org.Mm.eg.db},
      "YEAST" = {require(org.Sc.sgd.db); org.Sc.sgd.db}
   )

   # featureData <- AnnotationDbi::select(genome_data, keys = annot.keys,
   #                                      columns = annot.cols,
   #                                      keytype = keytype)

   featureData <- data.table() %>% .[, (keytype) := annot.keys]
   for (icol in annot.cols) {
      featureData <- AnnotationDbi::mapIds(genome_data,
                                           keys = annot.keys,
                                           column = icol,
                                           keytype = keytype,
                                           multiVals = AnnDbi.multiVals) %>%
         lapply(., function(x) paste(x, collapse = ";!") ) %>%
         do.call(rbind, .) %>%
         data.table(., keep.rownames = keytype) %>%
         setnames(., c(keytype, icol)) %>%
         merge(featureData, ., by = keytype)
   }


   featureData %<>% .[order(get(annot.cols[1]), decreasing = FALSE), ]
   ## use match !
   sel <- match(annot.keys, featureData[, ..keytype] %>% unlist)
   featureD <- featureData[sel, ]
   return(featureD)
}
