# reassemble_fig7.R — Figure 7 with new panel order: (a diagnosis | b BARD | c F3) / (d PLS | e secretome)
suppressMessages({library(ggplot2); library(patchwork); library(png); library(grid)})
O7 <- "/workspace/fig7_fig"
fg <- function(p, dir) wrap_elements(full = rasterGrob(readPNG(file.path(dir, p)), interpolate = TRUE))
p_diag <- fg("fig7c_diagnosis.png", O7)   # new a
p_bard <- fg("fig7d_bard.png", O7)        # new b
p_f3   <- fg("fig7e_f3_stratum.png", O7)  # new c
p_pls  <- fg("fig6g_pls_bridge.png", "/workspace")  # new d
p_sec  <- fg("fig6h_secretome.png", "/workspace")   # new e
top <- p_diag | p_bard | p_f3
bot <- (p_pls | p_sec) + plot_layout(widths = c(2, 3))
fig7 <- top / bot +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 14, family = "Liberation Sans"))
ggsave(file.path(O7, "fig7_main_neworder.svg"), fig7, width = 13, height = 8.6, bg = "white")
ggsave(file.path(O7, "fig7_main_neworder.png"), fig7, width = 13, height = 8.6, dpi = 300, bg = "white")
cat("FIG7_NEWORDER_DONE\n")
