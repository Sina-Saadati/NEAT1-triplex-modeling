source("analysis/00_setup.R")
db_unc <- "(((......)))..........((((((((.(((...))).)))))))).((((((((((..(((.(.(((((.....))).)).).))).))))))))))"  # RNAfold MFE, -35.00 kcal/mol
db_con <- ".....((((((((...........((((((.(((...))).))))))((((..........))))...))))))))........................."  # both regions single-stranded, -24.80 kcal/mol
stopifnot(nchar(db_unc) == 101, nchar(db_con) == 101)
db_pairs <- function(db) {
  ch <- strsplit(db, "")[[1]]; st <- integer(0); out <- list()
  for (i in seq_along(ch)) {
    if (ch[i] == "(") st <- c(st, i)
    else if (ch[i] == ")") { out[[length(out) + 1]] <- c(st[length(st)], i); st <- st[-length(st)] }
  }
  data.frame(i = sapply(out, `[`, 1), j = sapply(out, `[`, 2))
}
read_ct <- function(path) {                       # RhoFold+ (no MSA) 2D output
  ln <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8")); ln <- ln[grepl("^[0-9]+\\s+[ACGUT]", ln)]
  d <- read.table(text = ln, header = FALSE, stringsAsFactors = FALSE); names(d) <- c("i", "base", "prev", "nxt", "j", "idx")
  d[d$j > d$i, c("i", "j")]
}
pairs_all <- list(db_pairs(db_unc), db_pairs(db_con), read_ct(find_file("ss")))
YMAX <- max(sapply(pairs_all, function(p) max((p$j - p$i) / 2))) + 3
cats <- c("Other", "U-rich loop ↔ other", "A-rich tract ↔ other", "U-rich loop ↔ A-rich tract")
classify <- function(p) {
  inU <- function(x) x >= U[1] & x <= U[2]; inA <- function(x) x >= A[1] & x <= A[2]
  p$cat <- ifelse((inU(p$i) & inA(p$j)) | (inA(p$i) & inU(p$j)), cats[4],
                  ifelse(inA(p$i) | inA(p$j), cats[3], ifelse(inU(p$i) | inU(p$j), cats[2], cats[1])))
  p$cat <- factor(p$cat, levels = cats); p[order(p$cat), ]
}
arc_df <- function(p) do.call(rbind, lapply(seq_len(nrow(p)), function(k) {
  r <- (p$j[k] - p$i[k]) / 2; t <- seq(0, pi, length.out = 60)
  data.frame(x = (p$i[k] + p$j[k]) / 2 + r * cos(t), y = r * sin(t), grp = k, cat = p$cat[k]) }))
arc_plot <- function(p, ttl, sub = NULL, legend = FALSE) {
  p <- classify(p); cat(ttl, "->", nrow(p), "pairs;", paste(names(table(p$cat)), table(p$cat), collapse = " | "), "\n")
  g <- ggplot(arc_df(p), aes(x, y, group = grp, color = cat, linewidth = cat)) + geom_path(lineend = "round") +
    annotate("rect", xmin = U[1] - .5, xmax = U[2] + .5, ymin = -7, ymax = -2, fill = col[["orange"]]) +
    annotate("rect", xmin = A[1] - .5, xmax = A[2] + .5, ymin = -7, ymax = -2, fill = col[["verm"]]) +
    annotate("segment", x = .5, xend = 101.5, y = 0, yend = 0, color = "black", linewidth = .5) +
    scale_color_manual(values = setNames(c("#BDBDBD", col[["orange"]], col[["blue"]], col[["verm"]]), cats), name = NULL, drop = FALSE) +
    scale_linewidth_manual(values = setNames(c(.7, 1.1, 1.1, 1.8), cats), guide = "none") +
    scale_x_continuous(breaks = c(1, 20, 40, 60, 80, 101)) +
    coord_cartesian(xlim = c(0, 102), ylim = c(-8, YMAX), expand = FALSE) + labs(title = ttl, subtitle = sub) + th +
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(),
          plot.title.position = "plot", plot.title = element_text(face = "bold", size = 17, hjust = 0),
          plot.subtitle = element_text(size = 14, hjust = 0, margin = margin(b = 4)),
          legend.text = element_text(size = 13), legend.position = "bottom")
  if (legend) g + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 1.6))) else g + guides(color = "none")
}
pA <- arc_plot(pairs_all[[1]], "Unconstrained MFE (−35.00 kcal/mol)")
pB <- arc_plot(pairs_all[[2]], "Constrained MFE (−24.80 kcal/mol)", sub = "Both regions single-stranded; ΔΔG = +10.2 kcal/mol (RNAfold estimate)")
pC <- arc_plot(pairs_all[[3]], "RhoFold+ (no MSA) predicted 2D structure", legend = TRUE)
fig <- (pA / pB / pC) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
save_fig(fig, "Fig_secondary_structure_arcs", 8.5, 7.5)