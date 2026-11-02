#!/usr/bin/env Rscript
# fig3S1 v2: complete QC supplement — QC metrics + doublet rates, pre/post-integration
# UMAPs, canonical marker dotplots (human + mouse)
suppressMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(tidyr) })
source("/workspace/fig1_theme.R")
OUT <- "/workspace/fig3_out"; FIG <- "/workspace/fig3_fig"

h <- readRDS(file.path(OUT, "human_sc_scored.rds"))
mo <- readRDS(file.path(OUT, "mouse_sc_scored.rds"))

## (a) QC + doublet rates
qc_h <- read.csv(file.path(OUT, "qc_human_sc_raw.csv")); qc_h$species <- "Human"
qc_m <- read.csv(file.path(OUT, "qc_mouse_sc_raw.csv")); qc_m$species <- "Mouse"
qc <- bind_rows(qc_h, qc_m)
dbl <- read.csv(file.path(OUT, "doublet_rates.csv"))
qc <- qc %>% left_join(dbl %>% select(sample, doublet_rate, cells_postqc), by = "sample")
qc_l <- qc %>% pivot_longer(c(median_nFeature, median_nCount, median_mt, doublet_rate),
                            names_to = "metric", values_to = "value")
s1a <- ggplot(qc_l, aes(x = sample, y = value, fill = species)) +
  geom_col(width = 0.7) + facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(Human = palette8[1], Mouse = palette8[2])) +
  theme_ipf_lf() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3.5),
                         legend.position = "top") +
  labs(x = NULL, y = NULL)

## (b) pre- vs post-integration UMAP (human shown by group; mouse by condition)
h <- RunUMAP(h, reduction = "pca", dims = 1:30, reduction.name = "umap_raw", verbose = FALSE)
mo <- RunUMAP(mo, reduction = "pca", dims = 1:30, reduction.name = "umap_raw", verbose = FALSE)
p_pre_h <- DimPlot(h, reduction = "umap_raw", group.by = "group", pt.size = 0.05,
                   cols = c(Healthy = "#C8CACB", MASLD = "#E2EA4D")) +
  theme_ipf_lf() + ggtitle("Human pre-integration") + theme(legend.position = "none")
p_post_h <- DimPlot(h, reduction = "umap", group.by = "group", pt.size = 0.05,
                    cols = c(Healthy = "#C8CACB", MASLD = "#E2EA4D")) +
  theme_ipf_lf() + ggtitle("Human post-Harmony") + theme(legend.position = "none")
grp_pal_m <- c(CCl4_peak = "#E2EA4D", CCl4_reg2wk = "#B1C4D7", CCl4_reg6wk = "#BFDA8E",
               `FAT-MASH_peak` = "#E2EA4D", `FAT-MASH_reg2wk` = "#B1C4D7", `FAT-MASH_reg6wk` = "#BFDA8E")
p_pre_m <- DimPlot(mo, reduction = "umap_raw", group.by = "condition", pt.size = 0.05, cols = grp_pal_m) +
  theme_ipf_lf() + ggtitle("Mouse pre-integration") + theme(legend.position = "none")
p_post_m <- DimPlot(mo, reduction = "umap", group.by = "condition", pt.size = 0.05, cols = grp_pal_m) +
  theme_ipf_lf() + ggtitle("Mouse post-Harmony") + theme(legend.position = "none")

## (c) marker dotplots
mk_h <- c("ALB", "KRT19", "PECAM1", "COL1A1", "MARCO", "OLR1", "TREM2", "LYZ", "S100A8",
          "ITGAX", "MS4A2", "CD3D", "CD79A", "MZB1", "PPBP", "HBA1")
mk_h <- mk_h[mk_h %in% rownames(h)]
ct_order_h <- c("Hepatocyte", "Cholangiocyte", "Endothelial/LSEC", "HSC/myofibroblast",
                "Kupffer", "SAMac (OLR1+)", "Monocyte/MF", "Neutrophil", "cDC", "Mast cell",
                "T/NK", "B cell", "Plasma", "Platelet", "Other")
h$celltype_f <- factor(h$celltype_label, levels = ct_order_h)
dp_h <- DotPlot(h, features = mk_h, group.by = "celltype_f") +
  scale_color_gradient(low = "#ECE9E2", high = "#B1C4D7") +
  theme_ipf_lf() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                         axis.text.y = element_text(size = 6)) +
  labs(x = NULL, y = NULL, title = "Human markers")
mk_m <- c("Alb", "Krt19", "Pecam1", "Col1a1", "Clec4f", "Trem2", "Cd9", "Lyz2", "S100a8",
          "Clec9a", "Siglech", "Ms4a2", "Cd3d", "Cd79a", "Hbb-bt")
mk_m <- mk_m[mk_m %in% rownames(mo)]
ct_order_m <- c("Hepatocyte", "Cholangiocyte", "Endothelial/LSEC", "HSC/myofibroblast",
                "Kupffer", "SAMac (Trem2+)", "Monocyte/MF", "Neutrophil", "cDC", "pDC",
                "Mast cell", "T/NK", "B cell", "Erythrocyte")
mo$celltype_f <- factor(mo$celltype_label, levels = ct_order_m)
dp_m <- DotPlot(mo, features = mk_m, group.by = "celltype_f") +
  scale_color_gradient(low = "#ECE9E2", high = "#B1C4D7") +
  theme_ipf_lf() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                         axis.text.y = element_text(size = 6)) +
  labs(x = NULL, y = NULL, title = "Mouse markers")

s1 <- (s1a) / ((p_pre_h | p_post_h) / (p_pre_m | p_post_m)) / (dp_h | dp_m) +
  plot_layout(heights = c(1.1, 1.6, 1))
save_panel(s1, file.path(FIG, "fig3S1_qc"), width = 12, height = 13)
cat("S1 v2 DONE\n")
