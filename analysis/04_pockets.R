source("analysis/00_setup.R")
pk <- bind_rows(mutate(parse_info(find_file("RhoFold_noMSA_clean_info")), model = "RhoFold+ (no MSA)"),
                mutate(parse_info(find_file("RhoFold_withMSA_clean_info")), model = "RhoFold+ (with MSA)"),
                mutate(parse_info(find_file("RNAComposer_clean_info")), model = "RNAComposer"),
                mutate(parse_info(find_file("MD_final_100ps_clean_info")), model = "MD frame (100 ps)"))
pk$lab <- paste0(pk$model, " – P", pk$pocket,
                 ifelse(pk$model == "MD frame (100 ps)" & pk$pocket == 2, " (A-rich tract)",
                        ifelse(pk$model == "RNAComposer" & pk$pocket == 1, " (Loop L4)", "")))
pk <- pk %>% arrange(drug) %>% mutate(lab = factor(lab, levels = lab)); print(as.data.frame(pk))
cols <- c("RhoFold+ (no MSA)" = col[["blue"]], "RhoFold+ (with MSA)" = col[["purple"]],
          "RNAComposer" = col[["green"]], "MD frame (100 ps)" = col[["verm"]])
n <- nrow(pk)
p <- ggplot(pk, aes(drug, lab, color = model)) +
  annotate("rect", xmin = 0.2, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "#E69F00", alpha = .10) +
  annotate("rect", xmin = 0.5, xmax = 0.72, ymin = -Inf, ymax = Inf, fill = "#D55E00", alpha = .12) +
  annotate("text", x = 0.35, y = n + 0.9, label = "difficult", family = FONT, size = 5) +
  annotate("text", x = 0.61, y = n + 0.9, label = "druggable", family = FONT, size = 5) +
  geom_segment(aes(x = 0, xend = drug, yend = lab), linewidth = 0.9) + geom_point(size = 4.5) +
  geom_text(aes(label = sprintf("%.3f", drug)), hjust = -0.45, family = FONT, size = 4.8, color = "black") +
  scale_color_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(0, 0.72), breaks = seq(0, 0.6, 0.2), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.6))) +
  labs(x = "Druggability score (0–1)", y = NULL) + th + theme(axis.text.y = element_text(size = 14))
save_fig(p, "Fig_pocket_ranked", 6.5, 5.4)