#!/usr/bin/env Rscript
# Human curation part 2: diagnose ambiguous clusters + myeloid subclustering + final labels
suppressMessages({ library(Seurat); library(dplyr); library(UCell) })
OUT <- "/workspace/fig3_out"
m <- readRDS(file.path(OUT, "human_sc_integrated.rds"))
m <- JoinLayers(m)

## (a) diagnostic gene expression for ambiguous clusters
diag_genes <- c("ALB", "APOA2", "COL1A1", "COL3A1", "ACTA2", "RGS5", "MYH11", "LUM", "DCN",
                "KRT19", "EPCAM", "PECAM1", "STAB2", "CLEC4G", "CD3D", "NKG7", "LYZ",
                "MARCO", "CD5L", "TREM2", "CD9", "SPP1", "OLR1", "S100A8", "FCGR3B",
                "MS4A2", "CPA3", "MZB1", "HBA1", "HBB", "RGS5", "CCL21", "PROX1", "FLT4")
diag_genes <- unique(diag_genes[diag_genes %in% rownames(m)])
amb <- c("14", "15", "19", "20")
sub <- subset(m, cells = rownames(m@meta.data)[m$seurat_clusters %in% amb])
avg <- AverageExpression(sub, features = diag_genes, group.by = "seurat_clusters",
                         layer = "data")$RNA
cat("== diagnostic avg expression (log-normalized) ==\n")
print(round(avg, 2))

## (b) myeloid subclustering: cl 7, 10, 11
mye <- subset(m, cells = rownames(m@meta.data)[m$seurat_clusters %in% c("7", "10", "11")])
cat("myeloid cells:", ncol(mye), "\n")
mye <- NormalizeData(mye, verbose = FALSE) |> FindVariableFeatures(nfeatures = 2000, verbose = FALSE)
mye <- ScaleData(mye, verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE)
suppressMessages(library(harmony))
mye <- RunHarmony(mye, group.by.vars = "sample", verbose = FALSE)
mye <- RunUMAP(mye, reduction = "harmony", dims = 1:20, verbose = FALSE)
mye <- FindNeighbors(mye, reduction = "harmony", dims = 1:20, verbose = FALSE) |>
       FindClusters(resolution = 0.4, verbose = FALSE)
mmarkers <- list(
  Kupffer = c("MARCO", "CD5L", "VSIG4", "TIMD4", "CLEC4F"),
  SAMac_OLR1 = c("OLR1", "TREM2", "CD9", "SPP1", "GPNMB", "LGALS3", "FABP5"),
  Monocyte_MF = c("LYZ", "S100A8", "S100A9", "CCR2", "PLAC8", "VCAN", "FCN1"),
  Neutrophil = c("CSF3R", "CXCR2", "FCGR3B", "MMP9", "CEACAM8"),
  cDC = c("CLEC9A", "BATF3", "XCR1", "CLEC10A", "CD1C"),
  pDC = c("SIGLECH", "CLEC4C", "LILRA4")
)
d <- GetAssayData(mye, layer = "data")
sc <- as.data.frame(ScoreSignatures_UCell(d, features = mmarkers, name = "_mk"))
mye@meta.data <- cbind(mye@meta.data, sc)
cl_avg <- sapply(paste0(names(mmarkers), "_mk"),
                 function(v) tapply(mye@meta.data[[v]], mye$seurat_clusters, mean, na.rm = TRUE))
zs <- scale(cl_avg)
cat("== myeloid subcluster z-scores ==\n"); print(round(zs, 2))
write.csv(round(zs, 3), file.path(OUT, "human_myeloid_subcluster_zscores.csv"))

Idents(mye) <- mye$seurat_clusters
mk <- FindAllMarkers(mye, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10, with_ties = FALSE)
for (cl in levels(factor(top$cluster)))
  cat("sub", cl, ":", paste(top$gene[top$cluster == cl], collapse = ", "), "\n")
write.csv(mk, file.path(OUT, "human_myeloid_subcluster_markers.csv"), row.names = FALSE)
saveRDS(mye, file.path(OUT, "human_myeloid.rds"))
cat("DONE\n")
