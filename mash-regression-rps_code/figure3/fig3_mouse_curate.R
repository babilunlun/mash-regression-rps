#!/usr/bin/env Rscript
# Mouse: apply curated cell-type map + myeloid subclustering to resolve SAMac
suppressMessages({ library(Seurat); library(dplyr); library(UCell) })
OUT <- "/workspace/fig3_out"
m <- readRDS(file.path(OUT, "mouse_sc_integrated.rds"))
m <- JoinLayers(m)

cur <- c("0" = "B_cell", "1" = "Endothelial_LSEC", "2" = "Endothelial_LSEC",
         "3" = "T_NK", "4" = "Hepatocyte", "5" = "T_NK", "6" = "T_NK",
         "7" = "Hepatocyte", "8" = "Neutrophil_granulocyte", "9" = "Dendritic",
         "10" = "B_cell", "11" = "Monocyte_MF", "12" = "T_NK", "13" = "Kupffer",
         "14" = "Erythrocyte", "15" = "T_NK", "16" = "T_NK", "17" = "Dendritic",
         "18" = "T_NK", "19" = "Endothelial_LSEC", "20" = "HSC_myofibroblast",
         "21" = "Cholangiocyte")
m$celltype <- unname(cur[as.character(m$seurat_clusters)])
# subtype detail
subt <- c("5" = "NK", "15" = "Treg", "16" = "NK", "18" = "gdT", "9" = "pDC")
m$subtype <- ifelse(as.character(m$seurat_clusters) %in% names(subt),
                    unname(subt[as.character(m$seurat_clusters)]), m$celltype)

## ---- myeloid subclustering (Kupffer/Mono/SAMac/DC/Neutrophil clusters) ----
mye <- subset(m, cells = rownames(m@meta.data)[m$seurat_clusters %in% c("8", "9", "11", "13", "17")])
cat("myeloid cells:", ncol(mye), "\n")
mye <- NormalizeData(mye, verbose = FALSE) |> FindVariableFeatures(nfeatures = 2000, verbose = FALSE)
mye <- ScaleData(mye, verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE)
suppressMessages(library(harmony))
mye <- RunHarmony(mye, group.by.vars = "sample", verbose = FALSE)
mye <- RunUMAP(mye, reduction = "harmony", dims = 1:20, verbose = FALSE)
mye <- FindNeighbors(mye, reduction = "harmony", dims = 1:20, verbose = FALSE) |>
       FindClusters(resolution = 0.4, verbose = FALSE)

mmarkers <- list(
  Kupffer = c("Marco", "Cd5l", "Vsig4", "Timd4", "Clec4f"),
  SAMac_Trem2 = c("Trem2", "Cd9", "Spp1", "Gpnmb", "Fabp5", "Lgals3", "Lpl"),
  Monocyte_MF = c("Ly6c2", "Plac8", "Chil3", "Ccr2", "Vcan", "Fn1"),
  Neutrophil = c("S100a8", "S100a9", "Csf3r", "Cxcr2", "Retnlg"),
  cDC = c("Clec9a", "Batf3", "Xcr1"),
  pDC = c("Siglech", "Bst2", "Cox6a2")
)
d <- GetAssayData(mye, layer = "data")
sc <- as.data.frame(ScoreSignatures_UCell(d, features = mmarkers, name = "_mk"))
mye@meta.data <- cbind(mye@meta.data, sc)
cl_avg <- sapply(paste0(names(mmarkers), "_mk"),
                 function(v) tapply(mye@meta.data[[v]], mye$seurat_clusters, mean, na.rm = TRUE))
zs <- scale(cl_avg)
print(round(zs, 2))
write.csv(round(zs, 3), file.path(OUT, "mouse_myeloid_subcluster_zscores.csv"))

# top markers per subcluster
Idents(mye) <- mye$seurat_clusters
mk <- FindAllMarkers(mye, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10, with_ties = FALSE)
for (cl in levels(factor(top$cluster)))
  cat("sub", cl, ":", paste(top$gene[top$cluster == cl], collapse = ", "), "\n")
write.csv(mk, file.path(OUT, "mouse_myeloid_subcluster_markers.csv"), row.names = FALSE)
saveRDS(mye, file.path(OUT, "mouse_myeloid.rds"))

# provisional assignment by z-argmax (reviewed against top markers after)
sub_map <- setNames(colnames(zs)[apply(zs, 1, which.max)], rownames(zs))
print(sub_map)
m$celltype_fine <- m$celltype
mye_ct <- unname(sub_map[as.character(mye$seurat_clusters)])
m$celltype_fine[colnames(mye)] <- mye_ct
print(table(m$celltype_fine, m$condition))
saveRDS(m, file.path(OUT, "mouse_sc_integrated.rds"))
write.csv(m@meta.data, file.path(OUT, "mouse_sc_metadata.csv"))
cat("DONE\n")
