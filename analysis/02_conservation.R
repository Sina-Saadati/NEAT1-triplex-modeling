source("analysis/00_setup.R")
seqs <- read_fasta_md(find_file("homologs_aligned"))
aln <- do.call(rbind, strsplit(unname(seqs), "")); rownames(aln) <- names(seqs)
stopifnot(ncol(aln) == 101, nrow(aln) == 19)
bases <- c("A", "C", "G", "U")
cnt <- sapply(bases, function(b) colSums(aln == b))
cons_pct <- 100 * apply(cnt, 1, max) / nrow(aln)          # gaps count as mismatches
cons_base <- bases[apply(cnt, 1, which.max)]
cat(sprintf("U-rich %.1f%% | A-rich %.1f%% | nt85-101 %.1f%% | all %.1f%%\n",
            mean(cons_pct[U[1]:U[2]]), mean(cons_pct[A[1]:A[2]]), mean(cons_pct[85:101]), mean(cons_pct)))

pos <- 1:101
reg <- ifelse(pos >= U[1] & pos <= U[2], "U-rich loop (nt 14–23)",
              ifelse(pos >= A[1] & pos <= A[2], "A-rich tract (nt 82–101)", "Other positions"))
df <- data.frame(pos, cons_pct, cons_base,
                 reg = factor(reg, levels = c("U-rich loop (nt 14–23)", "A-rich tract (nt 82–101)", "Other positions")))
fig <- ggplot(df, aes(pos, cons_pct, fill = reg)) + geom_col(width = 0.85) +
  geom_text(data = df[df$pos == 15, ], aes(label = paste0(cons_base, pos), y = cons_pct + 4), hjust = 1.1, vjust = 0,
            family = FONT, size = 5.5, show.legend = FALSE) +
  geom_text(data = df[df$pos == 19, ], aes(label = paste0(cons_base, pos), y = cons_pct + 4), hjust = -0.1, vjust = 0,
            family = FONT, size = 5.5, show.legend = FALSE) +
  annotate("text", x = mean(c(85, 101)), y = 105, label = "nt 85–101: invariant (100%)", family = FONT, size = 5.5, vjust = 0) +
  scale_fill_manual(values = c("U-rich loop (nt 14–23)" = col[["orange"]], "A-rich tract (nt 82–101)" = col[["verm"]],
                               "Other positions" = "#B8B8B8")) +
  scale_x_continuous(limits = c(0, 104), breaks = c(1, 20, 40, 60, 80, 101), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + labs(x = "Alignment position (nt)", y = "Consensus frequency (%)") +
  th + theme(legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 15),
             axis.text = element_text(size = 16), plot.margin = margin(10, 32, 8, 8))
save_fig(fig, "Fig_conservation", 10.5, 4.6)