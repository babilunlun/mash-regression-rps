#!/usr/bin/env Rscript
# Figure 3 panel e: spatial mapping — mouse Visium (RCTD weights + component scores)
# + human GeoMx ROI scores
suppressMessages({ library(dplyr); library(ggplot2); library(patchwork); library(tidyr) })
source("/workspace/fig1_theme.R")
OUT <- "/workspace/fig3_out"; FIG <- "/workspace/fig3_fig"
DATA <- "/workspace/fig3_data/GSE306327"

read_pos <- function(posdir) {
  pos <- read.csv(file.path(DATA, posdir, "spatial", "tissue_positions.csv"))
  data.frame(barcode = pos$barcode, x = pos$pxl_col_in_fullres, y = pos$pxl_row_in_fullres,
             tissue = pos$in_tissue)
}

plots <- list()
for (sn in c("peak", "regression")) {
  w <- read.csv(file.path(OUT, paste0("rctd_weights_", sn, ".csv")), row.names = 1, check.names = FALSE)
  sc <- read.csv(file.path(OUT, paste0("spot_scores_", sn, ".csv")))
  pos <- read_pos(sn)
  df <- pos %>% inner_join(sc, by = "barcode") %>%
    inner_join(w %>% tibble::rownames_to_column("barcode"), by = "barcode") %>%
    filter(tissue == 1)
  df$y <- -df$y  # flip for anatomical orientation
  for (v in c("reversal_prog", "SAMac (Trem2+)", "Monocyte_MF", "HSC_myofibroblast")) {
    if (!v %in% colnames(df)) next
    plots[[paste(sn, v)]] <- ggplot(df, aes(x = x, y = y, color = .data[[v]])) +
      geom_point(size = 1.6) +
      scale_color_gradientn(colors = c("#F5F5F5", ifelse(grepl("prog", v), "#B1C4D7", "#C7B0C7"),
                                       ifelse(grepl("prog", v), "#5B7FA6", "#7A5F8A"))) +
      theme_ipf_lf() + theme(axis.text = element_blank(), axis.title = element_blank(),
                             axis.ticks = element_blank(), axis.line = element_blank(),
                             panel.grid = element_blank(),
                             legend.title = element_text(size = 6), legend.text = element_text(size = 5),
                             legend.key.size = unit(0.3, "cm")) +
      ggtitle(paste(sn, "-", gsub("_prog", " score", v))) +
      labs(color = NULL) + coord_fixed()
  }
}
pe_mouse <- (plots[["peak reversal_prog"]] | plots[["regression reversal_prog"]]) /
            (plots[["peak SAMac (Trem2+)"]] | plots[["regression SAMac (Trem2+)"]]) /
            (plots[["peak HSC_myofibroblast"]] | plots[["regression HSC_myofibroblast"]])
save_panel(pe_mouse, file.path(FIG, "fig3e_mouse_visium"), width = 8, height = 9)

## human GeoMx
score <- read.csv(file.path(OUT, "geomx_roi_scores.csv"))
score$group <- sub("genotype: ", "", score$group)
gs <- score %>% pivot_longer(c(reversal, emergent), names_to = "comp", values_to = "score")
pe1 <- ggplot(gs, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
  facet_wrap(~comp) +
  scale_fill_manual(values = c(Healthy = "#C8CACB", MASH = "#E2EA4D")) +
  theme_ipf_lf() + labs(x = NULL, y = "ROI score (mean z)", fill = NULL) +
  ggtitle("Human GeoMx CK8+ ROIs (GSE287826)")
pe2 <- ggplot(score, aes(x = samac, y = reversal, color = group)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
  scale_color_manual(values = c(Healthy = "#C8CACB", MASH = "#E2EA4D")) +
  theme_ipf_lf() + labs(x = "SAMac marker score", y = "Reversal score", color = NULL) +
  ggtitle("Reversal-SAMac spatial covariation")
pe_human <- pe1 | pe2
save_panel(pe_human, file.path(FIG, "fig3e_human_geomx"), width = 9, height = 4)
cat("PANEL E DONE\n")
