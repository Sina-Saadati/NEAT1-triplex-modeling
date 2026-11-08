source("analysis/00_setup.R"); library(gridExtra); library(grid)
steps <- data.frame(n = 1:8, row = c(1,1,1,1,2,2,2,2), colp = c(1,2,3,4,1,2,3,4),
                    title = c("Target domain", "2D structure", "Conservation", "3D modeling", "Model comparison", "MD simulation", "Pocket detection", "Docking"),
                    body = c("NEAT1 ENE+A-rich\nnt 22643–22743\n(101 nt, GQ859162.1)", "ViennaRNA MFE\nfree vs constrained\n(ΔΔG estimate)",
                             "BLASTn + MUSCLE\n19 mammalian\nspecies", "RhoFold+ (±MSA)\ntrRosettaRNA\nRNAComposer",
                             "C1′ RMSD matrix\npLDDT\nbarnaba base pairs", "GROMACS 2021.3\nAMBER χOL3, TIP3P\n100 ps + 500 ps ext.",
                             "fpocketR\n(fpocket 4.0.3)\nstatic vs MD frame", "AutoDock Vina 1.2.7\nsingle-frame\nMM-GBSA"),
                    grp = c("a","a","a","b","b","c","d","d"))
steps$x <- c(1.2, 3.55, 5.9, 8.25)[steps$colp]; steps$y <- c(3.3, 0.95)[steps$row]
fills <- c(a = "#FCE9C2", b = "#CFE3F2", c = "#CDEBDD", d = "#F5D5C3"); w <- 1.95; h <- 1.75
arr <- arrow(length = unit(0.18, "cm"), type = "closed")
ln <- function(x, xe, y, ye, a = FALSE) geom_segment(aes(x = x, xend = xe, y = y, yend = ye), linewidth = 0.7, color = "#444444",
                                                     arrow = if (a) arr else NULL, inherit.aes = FALSE, data = data.frame(x, xe, y, ye))
wf <- ggplot(steps) +
  geom_rect(aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = grp), color = "#444444", linewidth = 0.6) +
  geom_text(aes(x = x, y = y + 0.52, label = paste0(n, ". ", title)), family = FONT, fontface = "bold", size = 4.8) +
  geom_text(aes(x = x, y = y - 0.22, label = body), family = FONT, size = 4.2, lineheight = 0.95) +
  scale_fill_manual(values = fills, guide = "none") +
  ln(c(2.175, 4.525, 6.875), c(2.575, 4.925, 7.275), 3.3, 3.3, TRUE) + ln(c(2.175, 4.525, 6.875), c(2.575, 4.925, 7.275), 0.95, 0.95, TRUE) +
  ln(8.25, 8.25, 2.425, 2.0) + ln(8.25, 1.2, 2.0, 2.0) + ln(1.2, 1.2, 2.0, 1.825, TRUE) +
  coord_cartesian(xlim = c(0, 9.4), ylim = c(0, 4.3), expand = FALSE) + theme_void()
ggsave(file.path(OUT_DIR, "Fig_workflow.png"), wf, device = agg_png, width = 8, height = 3.8, units = "in", dpi = 600)

tb <- data.frame(Model = c("RhoFold+ (no MSA)", "RhoFold+ (with MSA)", "trRosettaRNA", "RNAComposer"),
                 Approach = c("Language-model DL", "DL + 19-species MSA", "Co-evolution DL", "Fragment assembly"),
                 `Mean pLDDT` = c("61.2", "68.7", "41.1", "n/a"), `Backbone long bonds` = c("1", "13", "n.d.", "n.d."),
                 Pockets = c("6", "1", "0", "1"), `Best druggability` = c("0.141", "0.050", "–", "0.555"), check.names = FALSE)
ttm <- ttheme_minimal(base_family = FONT, base_size = 15,
                      core = list(fg_params = list(hjust = 0, x = 0.04), bg_params = list(fill = c("white", "#F3F6FA"))),
                      colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.04), bg_params = list(fill = "#DCE6F1")))
tg <- tableGrob(tb, rows = NULL, theme = ttm)
agg_png(tempfile(fileext = ".png"), width = 14, height = 6, units = "in", res = 100)
w_in <- convertWidth(sum(tg$widths), "in", valueOnly = TRUE) + 0.4; h_in <- convertHeight(sum(tg$heights), "in", valueOnly = TRUE) + 0.4
invisible(dev.off())
agg_png(file.path(OUT_DIR, "Table_models.png"), width = w_in, height = h_in, units = "in", res = 600)
grid.newpage(); grid.draw(tg); invisible(dev.off())