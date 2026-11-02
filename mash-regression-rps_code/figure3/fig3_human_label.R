#!/usr/bin/env Rscript
# Human: apply curated labels (whole-object clusters + myeloid subcluster remap)
suppressMessages({ library(Seurat); library(dplyr) })
OUT <- "/workspace/fig3_out"
h <- readRDS(file.path(OUT, "human_sc_integrated.rds"))
h <- JoinLayers(h)
mye <- readRDS(file.path(OUT, "human_myeloid.rds"))

# verify ambiguous clusters with diagnostic genes
diag_genes <- c("ALB", "COL1A1", "COL3A1", "ACTA2", "PECAM1", "STAB2", "HBA1", "HBB", "KRT19", "CD3D", "LYZ")
diag_genes <- diag_genes[diag_genes %in% rownames(h)]
sub <- subset(h, cells = rownames(h@meta.data)[h$seurat_clusters %in% c("14", "15", "19", "20")])
avg <- AverageExpression(sub, features = diag_genes, group.by = "seurat_clusters", layer = "data")$RNA
cat("== diagnostics cl14/15/19/20 ==\n"); print(round(avg, 2))

cur <- c("0" = "Hepatocyte", "1" = "Hepatocyte", "2" = "Hepatocyte", "3" = "Endothelial/LSEC",
         "4" = "Hepatocyte", "5" = "T/NK", "6" = "HSC/myofibroblast", "7" = "Myeloid",
         "8" = "Hepatocyte", "9" = "T/NK", "10" = "Myeloid", "11" = "Myeloid",
         "12" = "T/NK", "13" = "Cholangiocyte", "14" = "HSC/myofibroblast",
         "15" = "Endothelial/LSEC", "16" = "Endothelial/LSEC", "17" = "B cell",
         "18" = "HSC/myofibroblast", "19" = "Other", "20" = "Endothelial/LSEC",
         "21" = "Plasma", "22" = "Mast cell")
# cl14: COL1A1/COL3A1/ACTA2 all high (ALB is ambient-RNA prone) -> HSC/myofibroblast
# cl19: mixed low-signal -> Other; cl15/20: PECAM1+STAB2 -> Endothelial/LSEC (already set)
h$celltype_label <- unname(cur[as.character(h$seurat_clusters)])

sub_map <- c("0" = "Other", "1" = "cDC", "2" = "Hepatocyte", "3" = "Monocyte/MF",
             "4" = "Monocyte/MF", "5" = "Kupffer", "6" = "Kupffer", "7" = "Neutrophil",
             "8" = "SAMac (OLR1+)", "9" = "Monocyte/MF", "10" = "T/NK",
             "11" = "HSC/myofibroblast", "12" = "cDC", "13" = "Platelet")
h$celltype_label[colnames(mye)] <- unname(sub_map[as.character(mye$seurat_clusters)])
print(table(h$celltype_label, h$group))
saveRDS(h, file.path(OUT, "human_sc_integrated_tmp.rds"))
file.rename(file.path(OUT, "human_sc_integrated_tmp.rds"), file.path(OUT, "human_sc_integrated.rds"))
write.csv(h@meta.data, file.path(OUT, "human_sc_metadata.csv"))
cat("DONE\n")
