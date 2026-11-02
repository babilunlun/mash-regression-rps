#!/usr/bin/env Rscript
# Resume annotation: ScoreSignatures_UCell (AddModuleScore_UCell has a layer-matching bug)
suppressMessages({ library(Seurat); library(dplyr); library(UCell) })
OUT <- "/workspace/fig3_out"
args <- commandArgs(trailingOnly = TRUE)
which <- args[1]  # "mouse" or "human"

if (which == "mouse") {
  markers <- list(
    Hepatocyte = c("Alb", "Apoa2", "Ttr", "Serpina1c"),
    HSC_myofibroblast = c("Col1a1", "Col3a1", "Acta2", "Pdgfrb", "Des", "Rgs5"),
    Endothelial_LSEC = c("Pecam1", "Kdr", "Stab2", "Clec4g", "Flt4"),
    Kupffer = c("Marco", "Cd5l", "Vsig4", "Timd4", "Clec4f"),
    SAMac_Trem2 = c("Trem2", "Cd9", "Spp1", "Gpnmb", "Fabp5", "Lgals3"),
    Monocyte_MF = c("Ly6c2", "Plac8", "Chil3", "Ccr2"),
    T_NK = c("Cd3d", "Cd3e", "Cd8a", "Nkg7"),
    B_cell = c("Cd79a", "Ms4a1", "Cd19"),
    Cholangiocyte = c("Krt19", "Krt7", "Epcam", "Sox9"),
    Dendritic = c("Itgax", "Clec9a", "Batf3", "Siglech")
  )
  rds <- file.path(OUT, "mouse_sc_integrated.rds")
} else {
  markers <- list(
    Hepatocyte = c("ALB", "APOA2", "TTR", "SERPINA1"),
    HSC_myofibroblast = c("COL1A1", "COL3A1", "ACTA2", "PDGFRB", "DES", "RGS5"),
    Endothelial_LSEC = c("PECAM1", "KDR", "STAB2", "CLEC4G", "FLT4"),
    Kupffer = c("MARCO", "CD5L", "VSIG4", "TIMD4", "CLEC4F"),
    SAMac_OLR1 = c("OLR1", "TREM2", "CD9", "SPP1", "GPNMB", "LGALS3"),
    Monocyte_MF = c("LYZ", "S100A8", "S100A9", "CCR2", "PLAC8"),
    T_NK = c("CD3D", "CD3E", "CD8A", "NKG7"),
    B_cell = c("CD79A", "MS4A1", "CD19"),
    Plasma = c("MZB1", "JCHAIN", "IGHG1"),
    Cholangiocyte = c("KRT19", "KRT7", "EPCAM", "SOX9"),
    Dendritic = c("ITGAX", "CLEC9A", "BATF3", "SIGLECH")
  )
  rds <- file.path(OUT, "human_sc_integrated.rds")
}

m <- readRDS(rds)
m <- JoinLayers(m)
d <- GetAssayData(m, layer = "data")
sc <- ScoreSignatures_UCell(d, features = markers, name = "_mark")
sc <- as.data.frame(sc)
stopifnot(identical(rownames(sc), colnames(m)))
m@meta.data <- cbind(m@meta.data, sc)
score_cols <- colnames(sc)
cl_avg <- sapply(score_cols, function(v) tapply(m@meta.data[[v]], m$seurat_clusters, mean, na.rm = TRUE))
ct_map <- setNames(colnames(cl_avg)[apply(cl_avg, 1, which.max)], rownames(cl_avg))
m$celltype <- unname(ct_map[as.character(m$seurat_clusters)])
grp <- if (which == "mouse") "condition" else "group"
print(table(m$celltype, m@meta.data[[grp]]))
saveRDS(m, rds)
write.csv(as.data.frame(cl_avg), file.path(OUT, paste0(which, "_cluster_marker_scores.csv")))
write.csv(m@meta.data, file.path(OUT, paste0(which, "_sc_metadata.csv")))
saveRDS(ct_map, file.path(OUT, paste0(which, "_celltype_map.rds")))
cat("ANNOTATE DONE:", which, "\n")
