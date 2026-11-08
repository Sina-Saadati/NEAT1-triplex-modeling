source("analysis/00_setup.R")
# ---- RMSD matrix ----
ln <- trimws(readLines(find_file("Step16_RNAComposer_comparison"), warn = FALSE, encoding = "UTF-8"))
ln <- ln[grepl("^[A-Za-z].* vs .*:\\s*[0-9.]+\\s*$", ln)]
m <- regmatches(ln, regexec("^(.+?) vs (.+?):\\s*([0-9.]+)\\s*$", ln, perl = TRUE))
prs <- data.frame(a = sapply(m, `[`, 2), b = sapply(m, `[`, 3), rmsd = as.numeric(sapply(m, `[`, 4)))
stopifnot(nrow(prs) == 6)
lab_map <- c("RhoFold(noMSA)" = "RhoFold+\n(no MSA)", "RhoFold(withMSA)" = "RhoFold+\n(with MSA)",
             "trRosettaRNA" = "trRosettaRNA", "RNAComposer" = "RNAComposer")
prs$a <- unname(lab_map[prs$a]); prs$b <- unname(lab_map[prs$b]); stopifnot(!anyNA(prs$a), !anyNA(prs$b))
lv <- unname(lab_map); ia <- match(prs$a, lv); ib <- match(prs$b, lv)
hm <- data.frame(col = factor(lv[pmin(ia, ib)], levels = lv[1:3]), row = factor(lv[pmax(ia, ib)], levels = rev(lv[2:4])), rmsd = prs$rmsd)
rng <- range(hm$rmsd)
p_hm <- ggplot(hm, aes(col, row, fill = rmsd)) + geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.2f", rmsd), color = rmsd > mean(rng)), family = FONT, size = 7, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_x_discrete(labels = function(x) ifelse(x == "trRosettaRNA", "trRosetta\nRNA", x)) +
  scale_fill_viridis_c(option = "D", direction = -1, name = "C1′ RMSD (Å)", limits = c(floor(rng[1]) - 1, ceiling(rng[2]) + 1),
                       guide = guide_colorbar(barheight = unit(4.5, "cm"), barwidth = unit(0.6, "cm"))) +
  coord_fixed() + th + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
                             axis.text = element_text(size = 14, lineheight = 0.9), legend.title = element_text(face = "bold", size = 15),
                             legend.text = element_text(size = 14), plot.margin = margin(4, 4, 4, 4))
save_fig(p_hm, "Fig_RMSD_matrix", 7.6, 4.6)

# ---- pLDDT (RhoFold+ values are stored in the PDB B-factor column) ----
read_plddt_csv <- function(path) { ln <- trimws(readLines(path, warn = FALSE)); ln <- ln[grepl("^[0-9]+,", ln)]
d <- read.csv(text = ln, header = FALSE); names(d) <- c("res", "plddt"); d }
read_pdb_b <- function(path) { ln <- readLines(path, warn = FALSE); ln <- ln[grepl("^ATOM", ln)]
d <- aggregate(as.numeric(substr(ln, 61, 66)), list(res = as.integer(substr(ln, 23, 26))), mean); names(d) <- c("res", "plddt"); d }
tr  <- read_plddt_csv(find_file("plddt")); tr$model <- "trRosettaRNA"
rf1 <- read_pdb_b(file.path(DATA_DIR, "NEAT1_ENE_A_RhoFold.pdb")); rf1$model <- "RhoFold+ (no MSA)"
rf2 <- read_pdb_b(file.path(DATA_DIR, "NEAT1_ENE_A_RhoFold_withMSA.pdb")); rf2$model <- "RhoFold+ (with MSA)"
pl <- bind_rows(rf1, rf2, tr); stopifnot(!anyNA(pl$plddt), all(table(pl$model) == 101))
pl$model <- factor(pl$model, levels = c("RhoFold+ (no MSA)", "RhoFold+ (with MSA)", "trRosettaRNA"))
print(pl %>% group_by(model) %>% summarise(all = mean(plddt), U_rich = mean(plddt[res >= U[1] & res <= U[2]]), A_rich = mean(plddt[res >= A[1] & res <= A[2]])))
p_pl <- ggplot(pl, aes(res, plddt, color = model)) +
  annotate("rect", xmin = U[1] - .5, xmax = U[2] + .5, ymin = 0, ymax = 100, fill = col[["orange"]], alpha = .25) +
  annotate("rect", xmin = A[1] - .5, xmax = A[2] + .5, ymin = 0, ymax = 100, fill = col[["verm"]], alpha = .20) +
  annotate("text", x = mean(U), y = 102, label = "U-rich loop", vjust = 0, family = FONT, size = 5.5) +
  annotate("text", x = mean(A), y = 102, label = "A-rich tract", vjust = 0, family = FONT, size = 5.5) +
  geom_hline(yintercept = c(50, 70), linetype = "dashed", color = "#9A9A9A") + geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("RhoFold+ (no MSA)" = col[["blue"]], "RhoFold+ (with MSA)" = col[["purple"]], "trRosettaRNA" = col[["green"]]), name = NULL) +
  scale_x_continuous(limits = c(0, 102), breaks = c(1, 20, 40, 60, 80, 101), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25), expand = c(0, 0)) + coord_cartesian(clip = "off") +
  labs(x = "Nucleotide (position)", y = "pLDDT (score)") + th +
  theme(legend.position = "top", legend.text = element_text(size = 15), plot.margin = margin(8, 24, 8, 8))
save_fig(p_pl, "Fig_pLDDT", 9, 4.4)

# ---- Docking (Vina) and MM-GBSA ----
ln <- readLines(find_file("Step13_docking_final_summary"), warn = FALSE, encoding = "UTF-8")
get_score <- function(nm) as.numeric(sub(paste0("^", nm, "\\s*:\\s*(-?[0-9.]+).*$"), "\\1", grep(paste0("^", nm, "\\s*:"), ln, value = TRUE)[1]))
sc <- data.frame(ligand = c("Kanamycin A", "Quercetin", "Caffeine"), score = c(get_score("Kanamycin A"), get_score("Quercetin"), get_score("Caffeine")))
stopifnot(!anyNA(sc$score)); sc$ligand <- factor(sc$ligand, levels = sc$ligand)
ha <- c("Kanamycin A" = 33, "Quercetin" = 22, "Caffeine" = 14)   # heavy atoms (C18H36N4O11, C15H10O7, C8H10N4O2)
sc$le <- -sc$score / ha[as.character(sc$ligand)]
exp_dG <- -8.60   # reported quercetin–MALAT1 value (Step13_docking_quercetin_results.txt)
pa <- ggplot(sc, aes(ligand, score, fill = ligand)) + geom_hline(yintercept = 0, linewidth = 0.5) + geom_col(width = 0.65) +
  geom_hline(yintercept = exp_dG, linetype = "dashed", color = col[["grey"]]) +
  geom_text(aes(y = 0.25, label = paste0(minus(sprintf("%.2f", score)), "\n(LE ", sprintf("%.2f", le), ")")), vjust = 0, family = FONT, size = 5.2, lineheight = 0.95) +
  annotate("text", x = 0.45, y = exp_dG - 0.55, hjust = 0, vjust = 1, family = FONT, size = 4.8, lineheight = 0.95, label = "dashed line: reported ΔG,\nquercetin–MALAT1 (\u22128.60)") +
  scale_fill_manual(values = c("Kanamycin A" = col[["blue"]], "Quercetin" = col[["orange"]], "Caffeine" = "#9A9A9A"), guide = "none") +
  scale_y_continuous(limits = c(-11.6, 3.2), breaks = seq(-10, 0, 2), labels = minus, expand = c(0, 0)) +
  labs(x = NULL, y = "Vina score (kcal/mol)") + th
mm <- paste(readLines(find_file("Step14_MMGBSA_results"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
num <- function(pat) as.numeric(regmatches(mm, regexec(pat, mm, perl = TRUE))[[1]][2])
total <- num("binding\\)\\s*=\\s*([+-]?[0-9.]+)")
vals <- c(num("VDWAALS=([+-]?[0-9.]+)"), num("EEL=([+-]?[0-9.]+)"), num("EGB=([+-]?[0-9.]+)"), num("ESURF=([+-]?[0-9.]+)"))
stopifnot(!anyNA(c(total, vals)), abs(sum(vals) - total) < 0.05)
end <- cumsum(vals); start <- c(0, head(end, -1))   # waterfall: each bar starts where the previous ended
wf <- data.frame(x = 1:5, ymin = c(pmin(start, end), min(0, total)), ymax = c(pmax(start, end), max(0, total)),
                 val = c(vals, total), type = c(ifelse(vals < 0, "fav", "unfav"), "total"))
wf$lab_y <- ifelse(wf$val >= 0, wf$ymax + 1.6, wf$ymin - 1.6)
conn <- data.frame(x = 1:4 + 0.4, xend = 2:5 - 0.4, y = end)
pb <- ggplot(wf) + geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_segment(data = conn, aes(x = x, xend = xend, y = y, yend = y), color = col[["grey"]], linewidth = 0.6) +
  geom_rect(aes(xmin = x - 0.4, xmax = x + 0.4, ymin = ymin, ymax = ymax, fill = type)) +
  geom_text(aes(x = x, y = lab_y, label = minus(sprintf("%+.2f", val))), family = FONT, size = 5.5) +
  scale_fill_manual(values = c(fav = col[["blue"]], unfav = col[["verm"]], total = "#444444"), guide = "none") +
  scale_x_continuous(breaks = 1:5, limits = c(0.4, 5.6), expand = c(0, 0),
                     labels = expression(Delta*E[vdW], Delta*E[ele], Delta*E[GB], Delta*E[surf], Delta*G[bind])) +
  scale_y_continuous(limits = c(-31, 6), breaks = seq(-30, 0, 10), labels = minus, expand = c(0, 0)) +
  labs(x = NULL, y = "Contribution (kcal/mol)") + th
fig_dock <- (pa | pb) + plot_layout(widths = c(1, 1.25)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 24, family = FONT))
save_fig(fig_dock, "Fig_docking_MMGBSA", 11.5, 4.8)