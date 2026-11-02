#!/usr/bin/env Rscript
# Combined Figure 3 from panel PNGs
suppressMessages({ library(cowplot); library(ggplot2); library(png); library(grid) })
source("/workspace/fig1_theme.R")
FIG <- "/workspace/fig3_fig"
P1 <- "/mnt/results/Figure3"  # phase-1 panels

img <- function(path) {
  p <- png::readPNG(path)
  ggdraw() + draw_image(p)
}

a <- img(file.path(P1, "fig3a_human_umap.png"))
b <- img(file.path(P1, "fig3b_ucell_scores.png"))
c <- img(file.path(P1, "fig3c_mouse_dynamics.png"))
d <- img(file.path(FIG, "fig3d_cellchat.png"))
e1 <- img(file.path(FIG, "fig3e_mouse_visium.png"))
e2 <- img(file.path(FIG, "fig3e_human_geomx.png"))
f <- img(file.path(FIG, "fig3f_misty_cosmx.png"))

row1 <- plot_grid(a, b, ncol = 2, rel_widths = c(1, 1.4), labels = c("a", "b"), label_size = 14)
row2 <- plot_grid(c, d, ncol = 2, rel_widths = c(1, 1.4), labels = c("c", "d"), label_size = 14)
row3 <- plot_grid(e1, NULL, e2, ncol = 1, rel_heights = c(1.7, 0.05, 0.75), labels = c("e", "", ""), label_size = 14)
row4 <- plot_grid(f, labels = "f", label_size = 14)
combined <- plot_grid(row1, row2, row3, row4, ncol = 1, rel_heights = c(1, 1, 2.1, 0.85))
save_panel(combined, "/workspace/fig3_fig/fig3_combined", width = 16, height = 24)
cat("COMBINED DONE\n")
