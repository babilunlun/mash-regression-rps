#!/usr/bin/env Rscript
# Human: marker z-score matrix + FindAllMarkers for curation (run after integration RDS exists)
suppressMessages({ library(Seurat); library(dplyr); library(UCell) })
OUT <- "/workspace/fig3_out"
m <- readRDS(file.path(OUT, "human_sc_integrated.rds"))
m <- JoinLayers(m)

markers <- list(
  Hepatocyte = c("ALB", "APOA2", "TTR", "SERPINA1"),
  HSC_myofibroblast = c("COL1A1", "COL3A1", "ACTA2", "PDGFRB", "RGS5"),
  Endothelial_LSEC = c("PECAM1", "KDR", "STAB2", "CLEC4G", "FLT4"),
  Kupffer = c("MARCO", "CD5L", "VSIG4", "TIMD4", "CLEC4F"),
  SAMac_OLR1 = c("OLR1", "TREM2", "CD9", "SPP1", "GPNMB", "LGALS3"),
  Monocyte_MF = c("LYZ", "S100A8", "S100A9", "CCR2", "PLAC8", "VCAN"),
  Neutrophil = c("S100A8", "S100A9", "CSF3R", "CXCR2", "FCGR3B"),
  T_NK = c("CD3D", "CD3E", "CD8A", "NKG7"),
  B_cell = c("CD79A", "MS4A1", "CD19"),
  Plasma = c("MZB1", "JCHAIN", "IGHG1"),
  Cholangiocyte = c("KRT19", "KRT7", "EPCAM", "SOX9"),
  Dendritic = c("ITGAX", "CLEC9A", "BATF3", "SIGLECH"),
  Mast = c("MS4A2", "CPA3", "TPSAB1", "KIT")
)
d <- GetAssayData(m, layer = "data")
sc <- as.data.frame(ScoreSignatures_UCell(d, features = markers, name = "_mk"))
stopifnot(identical(rownames(sc), colnames(m)))
m@meta.data <- cbind(m@meta.data, sc)
cl_avg <- sapply(paste0(names(markers), "_mk"),
                 function(v) tapply(m@meta.data[[v]], m$seurat_clusters, mean, na.rm = TRUE))
zs <- scale(cl_avg)
write.csv(round(zs, 3), file.path(OUT, "human_cluster_marker_zscores.csv"))
for (i in seq_len(nrow(zs)))
  cat(sprintf("cl%s | %s\n", rownames(zs)[i], paste(sprintf("%s=%.1f", colnames(zs), zs[i, ]), collapse = " ")))

Idents(m) <- m$seurat_clusters
mk <- FindAllMarkers(m, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10, with_ties = FALSE)
for (cl in levels(factor(top$cluster)))
  cat("cl", cl, ":", paste(top$gene[top$cluster == cl], collapse = ", "), "\n")
write.csv(mk, file.path(OUT, "human_cluster_markers.csv"), row.names = FALSE)
saveRDS(m, file.path(OUT, "human_sc_integrated.rds"))
cat("DONE\n")
