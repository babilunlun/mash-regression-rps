#!/usr/bin/env Rscript
# fig6s6_secreted_subset_fig.R — Fig. S8 (local file Fig6-S6): secreted-subset RPS feasibility
# a: Δsubset vs Δfibrosis (GSE193066)  b: null dist, Δ~Δfibrosis  c: null dist, subset~NAS
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(svglite)})
# project theme (fig1_theme.R no longer on disk; reconstructed from the shared definition)
theme_ipf_lf <- function(base = 12) {
  theme_classic(base_size = base) +
    theme(text = element_text(family = "Liberation Sans"),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          axis.text   = element_text(color = "black", size = base - 1),
          axis.title  = element_text(color = "black", size = base, face = "bold"),
          axis.line   = element_line(color = "black", linewidth = 0.5),
          axis.ticks  = element_line(color = "black", linewidth = 0.4),
          plot.title    = element_text(face = "bold", size = base + 1, color = "black"),
          plot.subtitle = element_text(size = base - 3, color = "#666666"),
          legend.background = element_rect(fill = alpha("white", 0.75), color = NA))
}
TB <- "/workspace/fig6/tables"
OUT <- "/workspace/fig6_fig"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

pd <- fread(file.path(TB, "rps_secreted_subset_paired.csv"))
nl <- fread(file.path(TB, "rps_secreted_subset_nulls.csv"))
pd[, class := factor(class, levels = c("regressor", "stable", "progressor"))]
cls3 <- c(regressor = "#BFDA8E", stable = "#C8CACB", progressor = "#DCA7EB")

pa <- ggplot(pd, aes(dfib, dsub, color = class)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6, show.legend = FALSE) +
  scale_color_manual(values = cls3, name = NULL) +
  labs(title = "a  Secreted-subset score vs histological course",
       subtitle = "Spearman rho = -0.29, P = 0.030 (partial | dNAS: -0.29, P = 0.029)",
       x = "Delta fibrosis stage (GSE193066)", y = "Delta secreted-subset score") +
  theme_ipf_lf(9) + theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
                          legend.text = element_text(size = 7),
                          legend.key.size = unit(0.35, "cm"),
                          legend.background = element_rect(fill = alpha("white", 0.75), color = NA))

null_panel <- function(v, obs, tag, xlab, pct_lab) {
  d <- data.table(rho = v)
  ggplot(d, aes(rho)) +
    geom_histogram(bins = 40, fill = "#C8CACB", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = median(v), linetype = 3, color = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = obs, color = "#FF9400", linewidth = 0.9) +
    labs(title = tag, subtitle = pct_lab, x = xlab, y = "Count") +
    theme_ipf_lf(9)
}
pb <- null_panel(nl$null_dfib, -0.285,
                 "b  Null: 1,000 random 54-gene down-arm subsets",
                 "Delta~Delta fibrosis Spearman rho (GSE193066)",
                 "observed (orange) at 68th percentile")
pc <- null_panel(nl$null_nas, -0.203,
                 "c  Null: subset score vs disease activity",
                 "Subset~NAS Spearman rho (GSE135251)",
                 "observed (orange) at 2nd percentile")

fig <- pa | pb | pc
svglite(file.path(OUT, "Fig6-S6_secreted_subset.svg"), width = 10.5, height = 3.3)
print(fig); dev.off()
png(file.path(OUT, "Fig6-S6_secreted_subset.png"), width = 10.5, height = 3.3,
    units = "in", res = 300, type = "cairo")
print(fig); dev.off()
cat("FIG_DONE\n")
