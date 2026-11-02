#!/usr/bin/env Rscript
# Human consolidated label + score (faithful to fig3_human_label.R + fig3_human_score.R,
# correct cell-name indexing; skips AUCell which is not needed for leave-out)
suppressMessages({ library(Seurat); library(dplyr); library(UCell) })
OUT <- "/workspace/fig3_out"
h <- readRDS(file.path(OUT, "human_sc_integrated.rds"))
h <- JoinLayers(h)
mye <- readRDS(file.path(OUT, "human_myeloid.rds"))

cur <- c("0" = "Hepatocyte", "1" = "Hepatocyte", "2" = "Hepatocyte", "3" = "Endothelial/LSEC",
         "4" = "Hepatocyte", "5" = "T/NK", "6" = "HSC/myofibroblast", "7" = "Myeloid",
         "8" = "Hepatocyte", "9" = "T/NK", "10" = "Myeloid", "11" = "Myeloid",
         "12" = "T/NK", "13" = "Cholangiocyte", "14" = "HSC/myofibroblast",
         "15" = "Endothelial/LSEC", "16" = "Endothelial/LSEC", "17" = "B cell",
         "18" = "HSC/myofibroblast", "19" = "Other", "20" = "Endothelial/LSEC",
         "21" = "Plasma", "22" = "Mast cell")
h$celltype_label <- unname(cur[as.character(h$seurat_clusters)])

sub_map <- c("0" = "Other", "1" = "cDC", "2" = "Hepatocyte", "3" = "Monocyte/MF",
             "4" = "Monocyte/MF", "5" = "Kupffer", "6" = "Kupffer", "7" = "Neutrophil",
             "8" = "SAMac (OLR1+)", "9" = "Monocyte/MF", "10" = "T/NK",
             "11" = "HSC/myofibroblast", "12" = "cDC", "13" = "Platelet")
mye_lab <- unname(sub_map[as.character(mye$seurat_clusters)])
fi <- match(colnames(mye), colnames(h))
h$celltype_label[fi] <- mye_lab
cat("=== celltype_label x group ===\n"); print(table(h$celltype_label, h$group))

rev_g <- intersect(readLines(file.path(OUT, "genes_reversal_human.txt")), rownames(h))
eme_g <- intersect(readLines(file.path(OUT, "genes_emergent_human.txt")), rownames(h))
cat("human genes: reversal", length(rev_g), "emergent", length(eme_g), "\n")
d <- GetAssayData(h, layer = "data")
sc <- as.data.frame(ScoreSignatures_UCell(d, features = list(reversal = rev_g, emergent = eme_g),
                                          name = "_prog"))
stopifnot(identical(rownames(sc), colnames(h)))
h@meta.data <- cbind(h@meta.data, sc)
saveRDS(h, file.path(OUT, "human_sc_scored.rds"))
write.csv(h@meta.data, file.path(OUT, "human_sc_metadata_rebuild.csv"))

## ---- anchors vs deployed ucell_summary_human.csv ----
md <- h@meta.data
cat("total:", nrow(md), "(expect 70424)\n")
sam <- md[md$celltype_label == "SAMac (OLR1+)", ]
print(table(sam$group))
cat("SAMac total:", nrow(sam), "(expect 195 = 52 Healthy + 143 MASLD)\n")
cat("SAMac reversal medians: Healthy", round(median(sam$reversal_prog[sam$group == "Healthy"]), 4),
    "MASLD", round(median(sam$reversal_prog[sam$group == "MASLD"]), 4),
    "(expect 0.0891 / 0.1026)\n")
rk <- md %>% group_by(celltype_label) %>%
  summarise(med_masld = median(reversal_prog[group == "MASLD"]), .groups = "drop") %>%
  arrange(desc(med_masld))
cat("=== MASLD reversal ranking (expect cDC #1, SAMac #2) ===\n"); print(as.data.frame(rk))
cat("DONE\n")
