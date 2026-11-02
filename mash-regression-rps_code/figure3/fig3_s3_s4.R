#!/usr/bin/env Rscript
# Figure 3 S3: CellChat full network circle plots (peak vs reg6wk, both models)
# Figure 3 S4: all cell-type RCTD weight maps on both Visium slides
suppressMessages({ library(CellChat); library(dplyr); library(ggplot2); library(patchwork) })
source("/workspace/fig1_theme.R")
OUT <- "/workspace/fig3_out"; FIG <- "/workspace/fig3_fig"
DATA <- "/workspace/fig3_data/GSE306327"

## ---- S3 ----
cc <- readRDS(file.path(OUT, "cellchat_objects.rds"))
for (mdl in c("CCl4", "FAT-MASH")) {
  for (tp in c("peak", "reg6wk")) {
    obj <- cc[[paste(mdl, tp, sep = "_")]]
    for (w in c("count", "weight")) {
      mat <- if (w == "count") obj@net$count else obj@net$weight
      png(file.path(FIG, paste0("fig3S3_network_", gsub("-", "", mdl), "_", tp, "_", w, ".png")),
          width = 1600, height = 1600, res = 200)
      par(mar = c(1, 1, 3, 1))
      netVisual_circle(mat, vertex.weight = as.numeric(table(obj@idents)),
                       weight.scale = TRUE, label.edge = FALSE,
                       title.name = paste(mdl, tp, "- interaction", w))
      dev.off()
      svg(file.path(FIG, paste0("fig3S3_network_", gsub("-", "", mdl), "_", tp, "_", w, ".svg")),
          width = 8, height = 8)
      par(mar = c(1, 1, 3, 1))
      netVisual_circle(mat, vertex.weight = as.numeric(table(obj@idents)),
                       weight.scale = TRUE, label.edge = FALSE,
                       title.name = paste(mdl, tp, "- interaction", w))
      dev.off()
    }
  }
}
cat("S3 DONE\n")

## ---- S4 ----
read_pos <- function(posdir) {
  pos <- read.csv(file.path(DATA, posdir, "spatial", "tissue_positions.csv"))
  data.frame(barcode = pos$barcode, x = pos$pxl_col_in_fullres, y = -pos$pxl_row_in_fullres,
             tissue = pos$in_tissue)
}
s4 <- list()
for (sn in c("peak", "regression")) {
  w <- read.csv(file.path(OUT, paste0("rctd_weights_", sn, ".csv")), row.names = 1, check.names = FALSE)
  pos <- read_pos(sn)
  df <- pos %>% inner_join(w %>% tibble::rownames_to_column("barcode"), by = "barcode") %>% filter(tissue == 1)
  for (ct in colnames(w)) {
    s4[[paste(sn, ct)]] <- ggplot(df, aes(x = x, y = y, color = .data[[ct]])) +
      geom_point(size = 0.9) +
      scale_color_gradientn(colors = c("#F5F5F5", "#C7B0C7", "#7A5F8A")) +
      theme_ipf_lf() + theme(axis.text = element_blank(), axis.title = element_blank(),
                             axis.ticks = element_blank(), axis.line = element_blank(),
                             panel.grid = element_blank(), legend.key.size = unit(0.25, "cm"),
                             legend.title = element_text(size = 5), legend.text = element_text(size = 4),
                             plot.title = element_text(size = 7)) +
      ggtitle(paste0(sn, " - ", ct)) + labs(color = NULL) + coord_fixed()
  }
}
for (sn in c("peak", "regression")) {
  ps <- s4[grepl(paste0("^", sn), names(s4))]
  combined <- wrap_plots(ps, ncol = 4)
  save_panel(combined, file.path(FIG, paste0("fig3S4_all_celltypes_", sn)), width = 14, height = 3.2 * ceiling(length(ps) / 4))
}
cat("S4 DONE\n")
