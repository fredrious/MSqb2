
########################### normalisation diagnostics
testNormMethod <- function(dt,
                           colexp = c("Replicate", "Condition"),
                           variable = "Abundance",
                           group = NULL,
                           wespal = NULL) {
  if (variable == "Intensity") dt$Intensity <- log2(dt$Intensity)
  dt[, mn := lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = "Feature", .SDcols = variable]
  dt[, sd := lapply(.SD, function(x) sd(x, na.rm = TRUE)), by = "Feature", .SDcols = variable]

  vsntest <- ggplot(unique(dt[, .(mn, sd)]), aes(x = rank(mn), y = sd)) +
    stat_binhex(aes(fill = ..count..), bins = 60) +
    scale_fill_viridis() +
    ggsci::scale_color_tron() +
    geom_quantile(aes(colour = as.factor(..quantile..)), method = "rqss", lambda = 150, size = 1, show.legend = TRUE, ) +
    theme_bw() +
    labs(x = "Rank(means)", y = "Standard deviation") +
    guides(col = guide_legend(title = "Quantile\nregression fit\n(lambda = 150)"))



  pairs <- combn(sort(as.character(unique(dt$Replicate))), 2)
  nScatter <- ncol(pairs)

  tmp <- summary(unlist(dt[, ..variable]))
  range <- c(floor(tmp["Min."]), ceiling(tmp["Max."]))

  gp.sct <- list()
  gp.ma <- list()
  count <- 1

  form <- "Feature + Condition ~ Replicate"
  dtw <- dcast(dt, form, value.var = variable)


  # require(ggpmisc)
  for (k in 1:nScatter) {
    gp.sct[[k]] <- ggscatter(dtw,
      x = pairs[1, k], y = pairs[2, k],
      color = "Condition", palette = "jco",
      add.params = list(size = 1.5, alpha = 0.1),
      add = "reg.line",
      conf.int = FALSE,
      alpha = 0.3,
      ylim = range,
      xlim = range,
      cor.coef.coord = range,
      legend = c(0.9, 0.2),
      ggtheme = theme_bw()
    ) +
      geom_abline(intercept = 0, slope = 1, col = "grey20", alpha = 0.6, lty = 2) +
      stat_cor(aes_string(color = "Condition"), label.x = 3, hjust = -1)



    dtw[, M := get(pairs[1, k]) - get(pairs[2, k])]
    dtw[, A := 0.5 * (get(pairs[1, k]) + get(pairs[2, k]))]

    gp.ma[[k]] <- ggscatter(dtw,
      x = "A", y = "M",
      color = "Condition", palette = "jco",
      add.params = list(size = 1.5, alpha = 0.1),
      add = "loess",
      alpha = 0.2,
      xlim = range,
      legend = c(0.9, 0.2),
      ggtheme = theme_bw()
    ) +
      geom_hline(yintercept = 0, col = "darkgreen", alpha = 0.6, lty = 2)
  }


  ggd <-
    ggdensity(dt,
      x = "Abundance", facet.by = "Condition",
      add = "mean", rug = FALSE,
      color = "Replicate",
      fill = "Condition",
      alpha = 0.4,
      palette = "jco",
      legend = "bottom",
      ggtheme = theme_bw()
    ) +
    color_palette(palette = "Greys")


  ggb <-
    ggboxplot(dt,
      y = "Abundance", facet.by = "Condition",
      color = "Condition",
      fill = "Condition",
      legend = "bottom",
      alpha = 0.4,
      palette = "jco",
      ggtheme = theme_bw()
    )
  # median_iqr: returns the median and the error limits defined by the interquartile range.
  ggb <- add_summary(ggb,
    group = "Replicate",
    fun = "median_iqr",
    size = 1.5,
    shape = 3,
    # color = "Replicate",
    position = position_dodge(0.7),
    error.plot = "pointrange"
  )
  ggpar(ggb, xlab = "")




  # m <- matrix(NA, ncol(dtest), ncol(dtest))
  # m[t(upper.tri(m, diag = F))] <- ((nScatter+ncol(dtest)+1) : (ncol(dtest)^2))
  # m <- t(m)
  # m[lower.tri(m, diag = F)] <- seq(nScatter)
  # diag(m) <- (nScatter+1) : (nScatter+ncol(dtest))

  # grid.arrange(grobs = c(gp.sct, gp.nrm, gp.ma) , layout_matrix = m)
  glist <- arrangeGrob(grobs = c(gp.sct, gp.nrm, gp.ma), layout_matrix = m)


  return(list(glist, vsntest))
}
###########################
