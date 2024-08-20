

msqb_config_viz <- function(build.para,
                            config.viz.file = NULL,
                            silent = FALSE,
                            #### Global plot parameters ----
                            manual.color.palette = NULL,
                            save.format = "pdf",
                            highlight = NULL, # Genes/Proteins to be manually added to heatmaps and volcano plots
                            ####  QC plots ----
                            ## Box and Density plots
                            BoxDistPlot.mag = 2,
                            # Dist.col = "BioRep",
                            # Dist.fill = "BioRep",
                            ## FractionDensQC
                            FractionQC.w = 15,
                            FractionQC.h = 4,
                            ## MissingPerFracChnl
                            MissingPerFracChnl.fill.alpha = 0.4,
                            ## extra
                            pl.width = "auto",
                            pl.height = "auto",
                            suffix = NULL,
                            #### PCA ----
                            pca.label = "auto",
                            pca.topNperc = 5,
                            pca.dims = c(1, 2),
                            pca.width = 7,
                            pca.height = 5,
                            pca.labelsize = 2,
                            pca.legend.position = "right",
                            #### Volcano plots
                            label = "Genes",
                            oob.xlim = c(-8, 8),
                            oob.ylim = c(0, 15),
                            fdr.hlines = c(0.1),
                            fdr.cutoff = 0.05,
                            logfc.cutoff = 1,
                            facet.text.size = 9,
                            col.plt = NULL,
                            x.ticks = "auto",
                            shape = NULL, # volcano plot points shape
                            hl.color = "orange",
                            hl.shape = 1,
                            hl.size = 3,
                            hl.txt.size = 3,
                            hl.fontface = "bold", # "plain", "bold", "italic", "bold.italic"
                            hl.force = 7,
                            hl.min.segment.length = 0.5,
                            hl.label.size = 0.1,
                            repel.force = 8,
                            repel.label.size = 0.2,
                            repel.size = 4,
                            top.by = "P.Value",
                            topN = 15,
                            legend.position = "right",
                            shorten.volcano.name = FALSE,
                            interactivePlot = TRUE,
                            #### Heatmaps
                            row.tag = "Genes",
                            cluster.rows = TRUE,
                            heatmap_legend_side = "right",
                            annotation_legend_side = "right",
                            HMname = "HeatMap",
                            body.color.palette = "viridis",
                            body.color.palette.scaled = c("#215b99", "grey95", "#ad2b0e"),
                            order.rows.by = NULL,
                            arrange.annot.bars.by.ref = NULL,
                            row.fontsize = 7,
                            col.fontsize = 8,
                            ttl.fontsize = 10,
                            column_names_rot = -45,
                            topNhm = 20,
                            topNhm.ind = 30,
                            SigHM.width = 8,
                            SigHM.height = 7,
                            SigHM.ind.width = 5,
                            SigHM.ind.height = 5,
                            #### ProfilePlot
                            ProfilePlot.facet = "auto",
                            ProfilePlot.width = 6,
                            ProfilePlot.height = 9,
                            #### VennDiagram
                            venn.input.subset = "all",
                            venn.palette = "Set3",
                            venn.output.plot = "both",
                            venn.fdr.cutoff = 0.1,
                            venn.set_size = 0, # 0 means no label
                            venn.edge_size = 3, # line thickness
                            venn.label = "count", # count or percentage or both
                            venn.label_geom = "label", # label or text
                            venn.label_alpha = 0.4,
                            venn.label_color = "darkred" ,#label color
                            venn.width = 9,
                            venn.height = 6,
                            upset.width = 6,
                            upset.height = 5,
                            ...) {

  ## import build parameters
  list2env(mget(c(
    "Scripts.path", "Tables.path", "analysis.name",
    "add.date.tag"),
  envir = as.environment(build.para)
  ),
  envir = environment()
  )


  ## list all args manually passed to the function
  viz.args <- as.list(match.call(expand.dots = TRUE))
  ## !!! this sometimes results in strange 'name' class arguments.
  ## A temporary work-around would be ...
  viz.args.d <- as.list(substitute(list(...)))[-1L]
  nm <- intersect(names(viz.args), names(viz.args.d))
  viz.args[nm] <- NULL
  viz.args <- append(viz.args, viz.args.d)
  rm(viz.args.d)
  ## this must be properly debugged in the future!
  ## !!!

  viz.args[which(names(viz.args) == "")] <- NULL
  viz.args[which(names(viz.args) %in% c("config.viz.file", "build.para"))] <- NULL

  ## importing config parameters from config file
  conf.viz.path <- MSqb2:::.importConfigfile(
    conf.fl = config.viz.file,
    Scripts.path = Scripts.path,
    analysis.name = analysis.name,
    add.date.tag = add.date.tag,
    whichConf = "viz",
    conf.args = viz.args
  )

  ## after importing config parameters from config file, force all
  ## parameters in ellipsis
  list2env(
    mget(names(viz.args), envir = as.environment(viz.args)),
    envir = environment()
  )





  # add parameters which will be modified to a second list
  viz.args2 <- list()
  ## ...
  ## here is to be completed similar to msqb_config
  ## ...




  ## remove items of viz.args already exist in viz.args2
  viz.args[names(viz.args) %in% names(viz.args2)] <- NULL
  ## now merge the two viz.args
  allargs <- c(viz.args2, viz.args)

  ## after parameter check, save changes to config parameters in config.viz.file
  if (length(allargs) > 0) {
    for (i in 1:length(allargs)) {
      assign(names(allargs)[i], allargs[[i]])

      txtin <- paste(names(allargs)[i], "=", paste(deparse(allargs[[i]]), collapse = "")) # deparse(allargs[[i]])
      code <- readLines(conf.viz.path) %>% gsub("<-", "=", .)
      nl <- which(lapply(strsplit(gsub(" ", "", code), split = "=", ), "[", 1) == names(allargs)[[i]])
      if (length(nl) > 1) {
        stop(
          "Multiple matches for ", names(allargs)[[i]],
          ". \nEdit the config file and rerun the function!"
        )
      } else if (length(nl) == 1) {
        txtout <- code[nl]
        code.edited <- gsub(txtout, txtin, code, fixed = TRUE)
        writeLines(code.edited, conf.viz.path)
        rm(code.edited)
      } else if (length(nl) == 0) {
        write(txtin, file = conf.viz.path, append = TRUE)
      }
    }
  }


  res <- svDialogs::dlg_message(
    paste0(
      "Do you want to view/edit the visualization config file? \n",
      "(If yes, save changes after editting and before ",
      "continuing with the visualisation step.)"
    ),
    "yesno"
  )$res
  if (res == "yes") file.edit(conf.viz.path)


  ## save configuration parameters in the Analyses/Tables directory.
  ## The file will be tagged by date and time for reproducibility purposes.
  paste0("_", format(Sys.time(), "%Y%m%d")) %>%
    paste0("VizParameters_", ., ".txt") %>%
    file.path(file.path(Tables.path, "Config_Parameters"), .) %>%
    file.copy(from = conf.viz.path, to = ., overwrite = FALSE)

  options(warn = -1)
  rm(
    allargs, viz.args2, viz.args,
    txtin, txtout, nl, code, i, res, config.viz.file
  )
  options(warn = 0)

  # some parameters might be of class expression or call. They must be evaluated!
  viz.conf.vars <- mget(ls(name = environment()), envir = environment()) %>%
    lapply(., eval)

  return(viz.conf.vars)
}
