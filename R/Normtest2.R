
########################### normalisation diagnostics
Normtest <- function(dt, colexp, var, group = NULL, wespal = NULL) {
  col.out <- setdiff(c("Abundance", "Intensity"), var)
  dt <- dt[, -..col.out]

  cols.fix <- names(dt[, !c(..colexp, ..var)]) # fixed columns
  form <- paste0(paste(cols.fix, collapse = "+"), "~", paste(colexp, collapse = "+")) # define casting formula
  dt.w <- dcast.data.table(dt, formula = form, value.var = eval(var)) # cast input data -- wide format
  cols.wid <- names(dt.w[, !..cols.fix]) # added columns after casting
  mat <- as.matrix(dt.w[, ..cols.wid]) # define matrix from expanded columns

  ## for a DIY meanSdPlot
  mu <- apply(mat, 1, function(x) mean(x, na.rm = TRUE))
  sd <- apply(mat, 1, function(x) sd(x, na.rm = TRUE))
  dtest <- data.table(cbind(mat, mu = rank(mu), sd))
  ## datatable solution!!! test this
  # dtest <- dt.w[,..cols.wid]
  # dtest[, `:=` (mu = rowMeans(.SD), sd = sd(.SD)), .SDcols = cols.wid]
  ##


  ## test vsn - QC
  vsntest <- ggplot(dtest, aes(x = mu, y = sd)) +
    stat_binhex(aes(fill = ..count..), bins = 60) +
    scale_fill_gradientn(colours = wes_palette(10, name = "Royal1", type = "continuous")) +
    geom_smooth(se = FALSE, method = lm, size = 1, col = "deepskyblue", alpha = 0.5) +
    geom_smooth(size = 1, col = "gold1", alpha = 0.5) +
    theme_bw() +
    labs(x = "Rank(means)", y = "Standard deviation")


  dtest <- as.data.table(mat)
  pairs <- combn(names(dtest), 2)
  nScatter <- ncol(pairs)

  tmp <- apply(dtest, 2, summary)
  range <- c(floor(min(tmp["Min.", ])), ceiling(max(tmp["Max.", ])))


  gp.sct <- list()
  gp.ma <- list()
  count <- 1

  for (k in 1:nScatter) {
    gp.sct[[k]] <-
      ggplot(dt.w, aes_string(x = pairs[1, k], y = pairs[2, k])) +
      geom_hex(aes(alpha = ..count.., fill = get(group)), bins = 50) +
      geom_abline(intercept = 0, slope = 1, col = "grey", alpha = 0.6) +
      geom_smooth(se = FALSE, method = lm, aes(col = get(group)), size = 1, alpha = 0.3) +
      stat_poly_eq(
        formula = y ~ x,
        aes(col = get(group), label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE
      ) +
      scale_fill_manual(values = c("#C29E0E", "#6D1314", "#7E8251"), name = "Condition") +
      scale_colour_manual(values = c("#C29E0E", "#6D1314", "#7E8251")) +
      scale_x_continuous(limits = range) +
      scale_y_continuous(limits = range) +
      guides(col = FALSE) +
      theme_bw() +
      theme(
        legend.position = c(0.9, 0.2),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "white", colour = "white")
      )



    MM <- dtest[, get(pairs[1, k])] - dtest[, get(pairs[2, k])]
    AA <- (dtest[, get(pairs[1, k])] + dtest[, get(pairs[2, k])]) / 2
    dma <- data.table("M" = MM, "A" = AA, dt.w[, ..group])

    gp.ma[[count]] <-
      ggplot(dma, aes(x = A, y = M)) +
      geom_hex(aes(fill = get(group), alpha = ..count..), bins = 40) +
      guides(alpha = FALSE) +
      geom_hline(yintercept = 0, col = "darkgreen", alpha = 0.6) +
      scale_fill_manual(values = c("#C29E0E", "#6D1314", "#7E8251"), name = "Condition") +
      scale_colour_manual(values = c("#C29E0E", "#6D1314", "#7E8251")) +
      theme_bw() +
      ggtitle(paste0(pairs[1, k], " v. ", pairs[2, k])) +
      scale_x_continuous(limits = range) +
      theme(
        legend.position = c(0.9, 0.2),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "white", colour = "white")
      )
    count <- count + 1
  }



  # man.palette = c("WT" = "#7E8251", "HET" = "#C29E0E", "KO" = "#6D1314")


  gp.nrm <- list()
  for (k in 1:ncol(dtest)) {
    dte <- data.table(dtest, dt.w[, ..group])

    gp.nrm[[k]] <-
      ggplot(dte, aes_string(colnames(dte)[k], color = eval(group))) +
      geom_histogram(position = "identity", binwidth = .6, aes(y = ..density.., fill = get(group)), alpha = 0.4) +
      geom_density() +
      guides(fill = FALSE) +
      scale_x_continuous(limits = range) +
      scale_fill_manual(values = c("#C29E0E", "#6D1314", "#7E8251"), name = "Condition") +
      scale_colour_manual(values = c("#C29E0E", "#6D1314", "#7E8251")) +
      theme_bw() +
      theme(
        legend.position = c(0.9, 0.2),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "white", colour = "white")
      )
  }



  m <- matrix(NA, ncol(dtest), ncol(dtest))
  m[t(upper.tri(m, diag = F))] <- ((nScatter + ncol(dtest) + 1):(ncol(dtest)^2))
  m <- t(m)
  m[lower.tri(m, diag = F)] <- seq(nScatter)
  diag(m) <- (nScatter + 1):(nScatter + ncol(dtest))

  glist <- arrangeGrob(grobs = c(gp.sct, gp.nrm, gp.ma), layout_matrix = m)

  return(list(glist, vsntest))
}



## display ggpubr paletes
## example 1: dis.cols(ggpubr::get_palette(palette = "futurama", 8),8)
## example 2: dis.cols(get_palette(palette =  c( "#FF0000", "#808080"), 6))
dis.cols <- function(a) image(1:length(a), 1, as.matrix(1:length(a)), col = a, xlab = "", ylab = "")

# for (plt in ggpubr:::.ggscipal()) {
#   for (j in 3:12) {
#     dis.cols(get_palette(palette =  plt, 6))
#
#   }
# }
