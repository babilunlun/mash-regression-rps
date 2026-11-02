#!/usr/bin/env Rscript
# Fig5-S5: cell-type attribution of the 3 candidate targets in mouse + human liver scRNA-seq.
# 2x2: (a) mouse expression dotplot, (b) human expression dotplot,
#      (c) mouse program-coupling heatmap, (d) human program-coupling heatmap.
.libPaths(c("/workspace/fig3_rlib/library", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})

TABS <- "/mnt/results/Figure5/tables"
OUT <- "/workspace/fig5s5"
dir.create(OUT, showWarnings = FALSE)

## ---------- load ----------
m_expr <- read.csv(file.path(TABS, "target3_mouse_expr_by_celltype.csv"), stringsAsFactors = FALSE)
h_expr <- read.csv(file.path(TABS, "target3_human_expr_by_celltype.csv"), stringsAsFactors = FALSE)
m_cor  <- read.csv(file.path(TABS, "target3_mouse_celllevel_cor.csv"), stringsAsFactors = FALSE)
h_cor  <- read.csv(file.path(TABS, "target3_human_celllevel_cor.csv"), stringsAsFactors = FALSE)
m_loo  <- read.csv(file.path(TABS, "target3_mouse_celllevel_cor_emergentLO.csv"), stringsAsFactors = FALSE)
h_loo  <- read.csv(file.path(TABS, "target3_human_celllevel_cor_emergentLO.csv"), stringsAsFactors = FALSE)

## ---------- expression dotplots ----------
gene_lab <- c("Trp53inp1" = "Trp53inp1", "Slamf6" = "Slamf6", "Fn3k" = "Fn3k",
              "TP53INP1" = "TP53INP1", "SLAMF6" = "SLAMF6", "FN3K" = "FN3K")

mk_dot <- function(d, title) {
  d <- d %>% mutate(celltype = reorder(celltype, -mean_expr))
  ggplot(d, aes(x = celltype, y = gene)) +
    geom_point(aes(size = pct_expr, color = mean_expr), alpha = 0.75) +
    scale_color_gradient(low = "grey90", high = "#0279EE", name = "mean\nexpr") +
    scale_size_continuous(range = c(0.2, 7), name = "%\nexpr") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(text = element_text(family = "Liberation Sans"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
          axis.text.y = element_text(face = "italic", size = 8),
          plot.title = element_text(size = 10, face = "bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7), legend.title = element_text(size = 7),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
}

pa <- mk_dot(m_expr, "a  Mouse liver scRNA-seq (n=27,202)")
pb <- mk_dot(h_expr, "b  Human liver scRNA-seq (n=70,424)")

## ---------- coupling heatmaps (LOO-corrected emergent) ----------
prep_cor <- function(cor, loo) {
  cor <- cor %>% filter(score != "emergent_prog")  # drop circular emergent
  loo <- loo %>% mutate(score = "emergent_prog")
  bind_rows(cor, loo) %>%
    mutate(score = factor(score, levels = c("reversal_prog", "emergent_prog",
                                            "SAMac_Trem2_mark", "SAMac_score", "HSC_myofibroblast_mark", "HSCact")),
           stratum = factor(stratum, levels = c("all", "myeloid", "SAMac", "HSC")))
}

score_lab <- c("reversal_prog" = "reversal", "emergent_prog" = "emergent (LOO)",
               "SAMac_Trem2_mark" = "SAMac", "SAMac_score" = "SAMac",
               "HSC_myofibroblast_mark" = "HSC act.", "HSCact" = "HSC act.")

mk_heat <- function(d, title) {
  d <- d %>% filter(!is.na(score)) %>%
    mutate(score = score_lab[as.character(score)],
           row = paste(gene, "|", score),
           star = ifelse(!is.na(p) & p < 0.001, "***",
                  ifelse(!is.na(p) & p < 0.01, "**",
                  ifelse(!is.na(p) & p < 0.05, "*", ""))))
  # order rows: gene then score
  score_ord <- c("reversal", "emergent (LOO)", "SAMac", "HSC act.")
  d <- d %>% mutate(score = factor(score, levels = score_ord)) %>%
    arrange(gene, score) %>% mutate(row = factor(row, levels = unique(row)))
  ggplot(d, aes(x = stratum, y = row)) +
    geom_tile(aes(fill = rho), color = "white", alpha = 0.9) +
    geom_text(aes(label = paste0(ifelse(is.na(rho), "NA", sprintf("%.2f", rho)), star)), size = 2.2) +
    scale_fill_gradient2(low = "#A7C0DF", mid = "white", high = "#FF9400",
                         midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish, name = "rho") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(text = element_text(family = "Liberation Sans"),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(size = 10, face = "bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7), legend.title = element_text(size = 7),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
}

pc <- mk_heat(prep_cor(m_cor, m_loo), "c  Mouse: cell-level coupling with programs")
pd <- mk_heat(prep_cor(h_cor, h_loo), "d  Human: cell-level coupling with programs")

fig <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1, 1.15))
ggsave(file.path(OUT, "Fig5-S5_cell_attribution.svg"), fig, width = 12, height = 9, device = svglite::svglite)
ggsave(file.path(OUT, "Fig5-S5_cell_attribution.png"), fig, width = 12, height = 9, dpi = 300)
cat("figure saved\n")
