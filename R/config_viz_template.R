config_viz_template <- function() {
## start config (don't touch this!)

  #### Global plot parameters ----
  manual.color.palette <- NULL # manual.color.palette = list(Condition = 'Set1', BioRep = 'Spectral')
  save.format <- "pdf"
  highlight <- NULL # Genes/Proteins to be manually added to heatmaps and volcano plots


  ####  QC plots ----
  ## Box and Density plots
  BoxDistPlot.mag <- 2
  # Dist.col = "BioRep"
  # Dist.fill = "BioRep"
  # Dist.fill.alpha = 0.6


  ## FractionDensQC
  FractionQC.w <- 15
  FractionQC.h <- 4


  ## MissingPerFracChnl
  MissingPerFracChnl.fill.alpha <- 0.4


  ## extra
  pl.width <- "auto"
  pl.height <- "auto"
  suffix <- NULL


  #### PCA
  pca.label <- "auto"
  pca.topNperc <- 5
  pca.dims <- c(1, 2)
  pca.width <- 7
  pca.height <- 5
  pca.labelsize <- 2
  pca.legend.position <- "right"


  #### Volcano plots
  label <- "Genes"
  fdr.hlines <- c(0.1)
  fdr.cutoff <- 0.05
  logfc.cutoff <- 1
  facet.text.size <- 9
  col.plt <- NULL
  x.ticks <- "auto"
  oob.xlim <-  c(-8, 8)
  oob.ylim <-  c(0, 15)
  shape <- NULL # volcano plot points shape
  hl.color <- "orange"
  hl.shape <- 1
  hl.size <- 3
  hl.txt.size <- 3
  hl.fontface <- "bold" # "plain"  "bold"  "italic"  "bold.italic"
  hl.force <- 7
  hl.min.segment.length <- 0.5
  hl.label.size <- 0.1
  repel.force <- 8
  repel.label.size <- 0.2
  repel.size <- 4
  top.by <- "P.Value"
  topN <- 10
  pval.sec.axis <- TRUE
  legend.position <- "right"
  shorten.volcano.name <- FALSE
  interactivePlot <- TRUE


  
  #### lfc bar plot
  lfc.width <- 12
  lfc.height <- 12



  #### Heatmaps
  row.tag <- "Genes"
  body.color.palette <- "viridis"
  scale.rows <- c(TRUE, FALSE)
  cluster.rows <- TRUE
  clustering_distance_rows <- "euclidean"
  clustering_method_rows <- "complete"
  clustering_distance_columns <- "euclidean"
  clustering_method_columns <- "complete"
  heatmap_legend_side <- "right"
  annotation_legend_side <- "right"
  HMname <- "HeatMap"
  order.rows.by <- NULL
  arrange.annot.bars.by.ref <- NULL
  row.fontsize <- 7
  col.fontsize <- 8
  ttl.fontsize <- 10
  column_names_rot <- -45
  topNhm <- 20
  topNhm.ind <- 30
  SigHM.width <- 8
  SigHM.height <- 7
  SigHM.ind.width <- 5
  SigHM.ind.height <- 5



  #### ProfilePlot
  ProfilePlot.facet <- "auto"
  topN.profileplot <- 20
  ProfilePlot.width <- 6
  ProfilePlot.height <- 9


  #### VennDiagram
  venn.input.subset <- "all"
  venn.palette <- "Set3"
  venn.variable <- "Protein"
  venn.output.plot <- "both"
  venn.fdr.cutoff <- 0.1
  venn.width <- 9
  venn.height <- 6
  upset.width <- 6
  upset.height <- 5

  #### extra parameters ----

## end config (don't touch this!)
}
