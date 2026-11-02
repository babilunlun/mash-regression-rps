#!/usr/bin/env Rscript
# make_fig7.R — Figure 7 (clinical benchmarking & translation) + Fig 6 re-layout (a-f)
# Fig 7: a PLS-NAFLD bridge | b secretome  /  c diagnosis | d BARD | e F3 stratum
# Fig 6: (a|b)/(c|d)/(e|f) — g/h moved to Fig 7
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(svglite);
                  library(png); library(grid)})
TH <- theme_classic(base_size = 10) +
  theme(text = element_text(family = "Liberation Sans"),
        axis.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "#666666"))
O7 <- "/workspace/fig7_fig"; dir.create(O7, showWarnings = FALSE, recursive = TRUE)
sv <- function(p, stem, w, h) {
  svglite(file.path(O7, paste0(stem, ".svg")), width = w, height = h); print(p); dev.off()
  png(file.path(O7, paste0(stem, ".png")), width = w, height = h, units = "in", res = 300,
      type = "cairo"); print(p); dev.off()
}

## ---------- panels c/d/e: benchmark trio (main-figure styling) ----------
d <- fread("/workspace/gse89632_rps_bard.csv")
d[, diag := factor(diag, levels = c("HC", "SS", "NASH"))]
pal <- c(HC = "#75A025", SS = "#FF9400", NASH = "#C44E52")
pc <- ggplot(d, aes(diag, RPS, fill = diag)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
  geom_jitter(width = 0.12, size = 0.9, alpha = 0.6) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(subtitle = "Kruskal P = 3.7e-5", x = NULL, y = "RPS") + TH
sv(pc, "fig7c_diagnosis", 3.4, 3.6)

pd_b <- d[!is.na(BARD)]
pd_b[, fib_f := factor(fib)]
pd <- ggplot(pd_b, aes(BARD, RPS)) +
  geom_boxplot(aes(group = BARD), outlier.shape = NA, alpha = 0.35, width = 0.55) +
  geom_jitter(aes(color = fib_f), width = 0.12, size = 1.1, alpha = 0.8) +
  scale_color_brewer(palette = "YlOrRd", name = "Fibrosis") +
  labs(subtitle = "Spearman rho = -0.37, P = 0.006", x = "BARD score", y = "RPS") + TH +
  theme(legend.position = "bottom", legend.text = element_text(size = 7),
        legend.title = element_text(size = 7.5))
sv(pd, "fig7d_bard", 3.6, 4.1)

f3 <- fread("/workspace/gse193066_f3_stratum.csv")
f3[, class := factor(class, levels = c("progressor", "stable", "regressor"))]
pe <- ggplot(f3, aes(class, RPS1, fill = class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
  geom_jitter(width = 0.12, size = 1.0, alpha = 0.65) +
  scale_fill_manual(values = c(progressor = "#C44E52", stable = "grey70",
                               regressor = "#75A025"), guide = "none") +
  labs(subtitle = "AUC = 0.754 (0.551-0.956), P = 0.042",
       x = NULL, y = "Baseline RPS (F3 stratum)") + TH +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8))
sv(pe, "fig7e_f3_stratum", 3.4, 3.6)

## ---------- assemble Figure 7: (a PLS | b secretome) / (c | d | e) ----------
fg <- function(p, dir) wrap_elements(full = rasterGrob(readPNG(file.path(dir, p)), interpolate = TRUE))
p_a <- fg("fig6g_pls_bridge.png", "/workspace")
p_b <- fg("fig6h_secretome.png", "/workspace")
p_c <- fg("fig7c_diagnosis.png", O7); p_d <- fg("fig7d_bard.png", O7); p_e <- fg("fig7e_f3_stratum.png", O7)
fig7 <- (p_a | p_b) / (p_c | p_d | p_e) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14, family = "Liberation Sans"))
ggsave(file.path(O7, "fig7_main.svg"), fig7, width = 13, height = 8.6, bg = "white")
ggsave(file.path(O7, "fig7_main.png"), fig7, width = 13, height = 8.6, dpi = 300, bg = "white")
cat("FIG7_DONE\n")

## ---------- re-assemble Fig 6 without g/h: (a|b)/(c|d)/(e|f) ----------
F6 <- "/mnt/results/Figure6/figures"
p_a <- fg("fig6a_schematic_v3.png", F6); p_b <- fg("fig6b_rps_course.png", F6)
p_c <- fg("fig6c_mouse_drug_response.png", F6); p_d <- fg("fig6d_stage_engagement.png", F6)
p_e <- fg("fig6e_roc.png", F6); p_f <- fg("fig6f_hcv_rygb_paired.png", F6)
fig6 <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14, family = "Liberation Sans"))
O6 <- "/workspace/fig6_fig"
ggsave(file.path(O6, "fig6_main.svg"), fig6, width = 13, height = 13.5, bg = "white")
ggsave(file.path(O6, "fig6_main.png"), fig6, width = 13, height = 13.5, dpi = 300, bg = "white")
cat("FIG6_DONE\n")
