#!/usr/bin/env Rscript
# Figure 3 panel d: CellChat peak vs reg6wk per model (CCl4, FAT-MASH)
suppressMessages({ library(Seurat); library(dplyr); library(CellChat) })
OUT <- "/workspace/fig3_out"
mo <- readRDS(file.path(OUT, "mouse_sc_scored.rds"))
mo <- JoinLayers(mo)

md <- mo@meta.data
md$cc_label <- md$celltype_label
md$cc_label[md$cc_label %in% c("cDC", "pDC")] <- "DC"
mo$cc_label <- md$cc_label

run_cc <- function(obj, cond) {
  sub <- subset(obj, cells = rownames(obj@meta.data)[obj@meta.data$condition == cond])
  keep <- names(which(table(sub$cc_label) >= 8))
  sub <- subset(sub, cells = rownames(sub@meta.data)[sub@meta.data$cc_label %in% keep])
  cat(cond, ": cells", ncol(sub), " groups:", paste(sort(unique(sub$cc_label)), collapse = ","), "\n")
  sub@meta.data$samples <- sub@meta.data$sample
  cc <- createCellChat(GetAssayData(sub, layer = "data"), meta = sub@meta.data, group.by = "cc_label")
  cc@DB <- CellChatDB.mouse
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 5)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

conds <- c("CCl4_peak", "CCl4_reg6wk", "FAT-MASH_peak", "FAT-MASH_reg6wk")
ccs <- lapply(conds, function(cd) run_cc(mo, cd))
names(ccs) <- conds
saveRDS(ccs, file.path(OUT, "cellchat_objects_tmp.rds"))
file.rename(file.path(OUT, "cellchat_objects_tmp.rds"), file.path(OUT, "cellchat_objects.rds"))

res <- list()
for (mdl in c("CCl4", "FAT-MASH")) {
  pk <- ccs[[paste0(mdl, "_peak")]]; rg <- ccs[[paste0(mdl, "_reg6wk")]]
  merged <- mergeCellChat(list(peak = pk, reg6wk = rg), add.names = c("peak", "reg6wk"))
  res[[mdl]] <- merged
  df_n <- compareInteractions(merged, show.legend = FALSE, measure = "count", return.data = TRUE)
  df_w <- compareInteractions(merged, show.legend = FALSE, measure = "weight", return.data = TRUE)
  write.csv(df_n, file.path(OUT, paste0("cellchat_diff_count_", mdl, ".csv")), row.names = FALSE)
  write.csv(df_w, file.path(OUT, paste0("cellchat_diff_weight_", mdl, ".csv")), row.names = FALSE)
  gg <- rankNet(merged, mode = "comparison", measure = "weight", sources.use = NULL,
                targets.use = NULL, do.stat = TRUE, return.data = TRUE)
  write.csv(gg$signaling.contribution, file.path(OUT, paste0("cellchat_ranknet_", mdl, ".csv")), row.names = FALSE)
}
saveRDS(res, file.path(OUT, "cellchat_merged_tmp.rds"))
file.rename(file.path(OUT, "cellchat_merged_tmp.rds"), file.path(OUT, "cellchat_merged.rds"))

mye <- c("SAMac (Trem2+)", "Monocyte/MF", "Kupffer", "Neutrophil", "DC")
for (mdl in c("CCl4", "FAT-MASH")) {
  for (tp in c("peak", "reg6wk")) {
    cc <- ccs[[paste0(mdl, "_", tp)]]
    df <- subsetCommunication(cc)
    df2 <- df[df$source %in% mye & df$target == "HSC/myofibroblast", ]
    write.csv(df2, file.path(OUT, paste0("cellchat_myeloid_to_hsc_", mdl, "_", tp, ".csv")), row.names = FALSE)
    cat(mdl, tp, "myeloid->HSC interactions:", nrow(df2), "\n")
  }
}
cat("DONE\n")
