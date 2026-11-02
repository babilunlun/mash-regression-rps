#!/usr/bin/env Rscript
# Fig3-S9 figure: human SAMac<->reversal leave-markers-out sensitivity (3 panels)
suppressMessages({ library(dplyr); library(ggplot2); library(tidyr); library(patchwork) })
source("/workspace/fig1_theme.R")
OUT <- "/workspace/fig3_out"

rankA <- read.csv(file.path(OUT, "fig3s9_celltype_ranking_leaveout_human.csv"))
corB  <- read.csv(file.path(OUT, "fig3s9_cor_full_vs_leaveout_human.csv"))
testD <- read.csv(file.path(OUT, "fig3s9_samac_masld_vs_healthy_human.csv"))

## ---- panel a: MASLD cell-type medians, full vs leave-out (SAMac highlighted) ----
pa_df <- rankA %>% filter(group == "MASLD") %>%
  mutate(celltype_label = reorder(celltype_label, median_lo)) %>%
  select(celltype_label, median_full, median_lo) %>%
  pivot_longer(-celltype_label, names_to = "geneset", values_to = "med") %>%
  mutate(geneset = recode(geneset, median_full = "Full (65 genes)",
                          median_lo = "Leave-out (61 genes)"))
pa <- ggplot(pa_df, aes(x = med, y = celltype_label, color = geneset)) +
  geom_line(aes(group = celltype_label), color = "grey75", linewidth = 0.5) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("Full (65 genes)" = "#DCA7EB",
                                "Leave-out (61 genes)" = "#0279EE")) +
  labs(x = "Median reversal score (MASLD)", y = NULL, color = NULL) +
  theme_ipf_lf() + theme(legend.position = "top", legend.text = element_text(size = 8.5))

## ---- panel b: per-group Spearman rho, full vs leave-out ----
pb_df <- corB %>%
  select(group, rho_full, rho_lo) %>%
  pivot_longer(-group, names_to = "geneset", values_to = "rho") %>%
  mutate(geneset = recode(geneset, rho_full = "Full", rho_lo = "Leave-out"),
         group = factor(group, levels = c("Healthy", "MASLD")))
pb <- ggplot(pb_df, aes(x = rho, y = group, color = geneset)) +
  geom_line(aes(group = group), color = "grey70", linewidth = 0.5) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c("Full" = "#DCA7EB", "Leave-out" = "#0279EE")) +
  labs(x = "Spearman rho (SAMac score vs reversal score, myeloid)", y = NULL, color = NULL) +
  theme_ipf_lf() + theme(legend.position = "top")

## ---- panel c: SAMac cluster Healthy vs MASLD, full vs leave-out ----
pc_df <- testD %>%
  select(full_healthy, full_masld, lo_healthy, lo_masld) %>%
  pivot_longer(everything(), names_to = c("geneset", "grp"), names_sep = "_", values_to = "med") %>%
  mutate(geneset = recode(geneset, full = "Full", lo = "Leave-out"),
         grp = factor(grp, levels = c("healthy", "masld")))
pc <- ggplot(pc_df, aes(x = grp, y = med, group = geneset, color = geneset)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
  scale_color_manual(values = c("Full" = "#DCA7EB", "Leave-out" = "#0279EE")) +
  scale_x_discrete(labels = c("healthy" = "Healthy", "masld" = "MASLD"), expand = c(0.3, 0.3)) +
  labs(x = NULL, y = "SAMac (OLR1+) median reversal score", color = NULL) +
  theme_ipf_lf() + theme(legend.position = "top")

fig <- (pa | pb | pc) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 13))
save_panel(fig, file.path(OUT, "fig3s9_samac_leaveout_human"), width = 12, height = 4.6)
cat("DONE\n")
