


# pEXPvLFC <-
#   ggplot(statLists[!is.na(adj.P.Val)], aes(x = AveExpr,  y = logFC, label1 = Gene, label2 = adj.P.Val, label3 = logFC)) +
#   geom_point(statLists[adj.P.Val > 0.05],  mapping = aes(alpha = abs(logFC), size = -log10(adj.P.Val)), color = "grey50") +
#   geom_point(statLists[adj.P.Val < 0.05],  mapping = aes(alpha = abs(logFC), size = -log10(adj.P.Val)), color = "red3") +
#   facet_grid(Comparison ~., scales = "free_x") +
#   geom_hline(yintercept = c(-1, 1), col = "grey20", alpha = 0.8, size = 0.3, linetype = "dashed") +
#   geom_hline(yintercept = 0, col = "grey20", alpha = 0.8, size = 0.3, linetype = "solid") +
#   theme_bw() +
#   labs(size = "Significance:\n-log10(FDR)", alpha = "Fffect size:\n|logFC|")
#
#
# pEXPvLFC.lbl <- pEXPvLFC +
#   geom_label_repel(data = statLists[adj.P.Val < 0.05] , col = "red3", aes(label = Gene ),
#                    force = 12,
#                    segment.size  = 0.2,
#                    # direction     = "x",
#                    # vjust = 2,
#                    show.legend = FALSE,
#                    label.size = 0.1,
#                    alpha = 0.7,
#                    point.padding = unit(0.5, "lines"),
#                    box.padding = unit(0.1, "lines"),
#                    segment.alpha = 0.8 ,
#                    segment.colour = "red3")
#
#
# pEXPvLFC.inta <- ggplotly(pEXPvLFC, tooltip = c("label1", "label2", "label3"))
#
#
# path <- file.path(adir, "Volcanos")
#
# htmlwidgets::saveWidget(as_widget(pEXPvLFC.inta), file = file.path(path, "pEXPvLFC_interactive.html"), selfcontained = TRUE)
#
# if (toupper(save.format) == "PNG") {
#   png(file.path(path, "EXPvLFC.png" ), width = 10, height = 8, units = "in", res = 150)
# }
# if (toupper(save.format) == "PDF") {
#   pdf(file.path(path, "EXPvLFC.pdf" ), width = 10, height = 8)
# }
# print(pEXPvLFC)
# dev.off()
#
#
# if (toupper(save.format) == "PNG") {
#   png(file.path(path, "EXPvLFC_labeled.png" ), width = 10, height = 8, units = "in", res = 150)
# }
# if (toupper(save.format) == "PDF") {
#   pdf(file.path(path, "EXPvLFC_labeled.pdf" ), width = 10, height = 8)
# }
# print(pEXPvLFC.lbl)
# dev.off()
