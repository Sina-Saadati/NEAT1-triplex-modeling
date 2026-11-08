source("analysis/00_setup.R")
dedupe <- function(df) df %>% distinct(t, .keep_all = TRUE) %>% arrange(t)   # restart-duplicated frames
smooth <- function(y, w = 11) { n <- length(y); h <- (w - 1) / 2
vapply(seq_len(n), function(i) mean(y[max(1, i - h):min(n, i + h)]), numeric(1)) }

rmsd <- dedupe(read_xvg(find_file("rmsd_500ps"), 2, c("t", "y")))
hb   <- dedupe(read_xvg(find_file("hbond_count_500ps"), 3, c("t", "y", "pairs")))
dist <- dedupe(read_xvg(find_file("distance_500ps"), 2, c("t", "y")))
cat(sprintf("Frames: %d | RMSD mean %.3f nm | H-bond mean %.2f (%d-%d) | COM dist %.3f +/- %.3f nm\n",
            nrow(rmsd), mean(rmsd$y), mean(hb$y), min(hb$y), max(hb$y), mean(dist$y), sd(dist$y)))

panel <- function(df, colr, ylab, steps = FALSE, hline = NULL) {
  df$s <- smooth(df$y)
  p <- ggplot(df, aes(t, y)) +
    { if (steps) geom_step(color = colr, alpha = .5, linewidth = .5) else geom_line(color = colr, alpha = .45, linewidth = .5) } +
    geom_line(aes(y = s), color = colr, linewidth = 1.3) +
    scale_x_continuous(expand = expansion(mult = .02)) + scale_y_continuous(expand = expansion(mult = .08)) +
    labs(x = NULL, y = ylab) + th
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = col[["grey"]], linewidth = .6)
  p
}
p1 <- panel(rmsd, col[["blue"]], "RMSD (nm)")
p2 <- panel(dist, col[["green"]], "COM distance (nm)")
p3 <- panel(hb, col[["verm"]], "H-bonds (number)", steps = TRUE, hline = mean(hb$y)) +
  scale_y_continuous(breaks = seq(2, 12, 2), expand = expansion(mult = .08)) + labs(x = "Extended-run time (ps)")
save_fig((p1 / p2 / p3) + plot_annotation(tag_levels = "a"), "Fig_MD_timeseries", 7.5, 9.5)

# RMSF (initial 100-ps run)
rf <- read_xvg(find_file("rmsf"), 2, c("nt", "rmsf"))
top <- rf %>% filter(nt >= A[1], nt <= A[2]) %>% slice_max(rmsf, n = 2) %>% arrange(nt)
ymax_rect <- max(rf$rmsf) + 0.045; y_label <- ymax_rect + 0.012
fig_rmsf <- ggplot(rf, aes(nt, rmsf)) +
  annotate("rect", xmin = U[1] - .5, xmax = U[2] + .5, ymin = -Inf, ymax = ymax_rect, fill = col[["orange"]], alpha = .30) +
  annotate("rect", xmin = A[1] - .5, xmax = A[2] + .5, ymin = -Inf, ymax = ymax_rect, fill = col[["verm"]], alpha = .22) +
  annotate("text", x = mean(U), y = y_label, label = "U-rich loop", vjust = 0, family = FONT, size = 5.5) +
  annotate("text", x = mean(A), y = y_label, label = "A-rich tract", vjust = 0, family = FONT, size = 5.5) +
  geom_hline(yintercept = mean(rf$rmsf), linetype = "dashed", color = col[["grey"]]) +
  geom_line(color = col[["blue"]], linewidth = 1.2) + geom_point(data = top, color = col[["verm"]], size = 3) +
  geom_text_repel(data = top, aes(label = paste("nt", nt)), family = FONT, size = 5.5, nudge_x = c(-7, 5),
                  nudge_y = 0.035, min.segment.length = 0, segment.color = col[["grey"]], seed = 1) +
  scale_x_continuous(expand = expansion(mult = .02)) + scale_y_continuous(expand = expansion(mult = c(.08, .32))) +
  labs(x = "Nucleotide (position)", y = "RMSF (nm)") + th
save_fig(fig_rmsf, "Fig_RMSF", 8, 4.2)