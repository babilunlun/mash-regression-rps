#!/usr/bin/env Rscript
# Figure 3 panel f (human part): CosMx LR spatial colocalization
# For key macrophage->stellate LR pairs: receptor expression in stellate cells
# adjacent to ligand+ macrophages vs distant stellate cells.
suppressMessages({ library(Seurat); library(dplyr); library(tidyr) })
OUT <- "/workspace/fig3_out"

lr_pairs <- list(
  `SPP1-CD44` = c("SPP1", "CD44"),
  `TGFB1-TGFBR1` = c("TGFB1", "TGFBR1"),
  `THBS1-CD47` = c("THBS1", "CD47"),
  `THBS1-CD36` = c("THBS1", "CD36"),
  `MIF-ACKR3` = c("MIF", "ACKR3"),
  `OSM-OSMR` = c("OSM", "OSMR"),
  `CCL2-CCR2` = c("CCL2", "CCR2"),
  `CXCL10-CXCR3` = c("CXCL10", "CXCR3")
)

res_all <- list()
for (slide in c("Slide1", "Slide2")) {
  f <- file.path("/workspace/fig3_data/GSE292268",
                 ifelse(slide == "Slide1", "GSM8855042_Slide1_SeuratObject.RDS", "GSM8855043_Slide2_SeuratObject.RDS"))
  s <- readRDS(f)
  ct_col <- grep("nbclust.*clusters$", colnames(s@meta.data), value = TRUE)[1]
  s$ct <- s@meta.data[[ct_col]]
  coords <- s@meta.data[, c("x_slide_mm", "y_slide_mm")]
  colnames(coords) <- c("x", "y")
  d <- GetAssayData(s, layer = "data")
  mac_types <- c("Inflammatory.macrophages", "Non.inflammatory.macrophages")
  ste <- colnames(s)[s$ct == "Stellate.cells"]
  mac <- colnames(s)[s$ct %in% mac_types]
  cat(slide, ": stellate", length(ste), " macrophages", length(mac), "\n")

  # distance matrix stellate x macrophage (mm -> um)
  dm <- as.matrix(dist(rbind(coords[ste, ], coords[mac, ]))) * 1000
  dm <- dm[seq_along(ste), length(ste) + seq_along(mac), drop = FALSE]

  for (lr in names(lr_pairs)) {
    lg <- lr_pairs[[lr]][1]; rc <- lr_pairs[[lr]][2]
    if (!lg %in% rownames(d) || !rc %in% rownames(d)) next
    lig_expr <- d[lg, mac]
    lig_pos <- mac[lig_expr > 0]
    if (length(lig_pos) < 10) next
    # stellate cells within 50 um of >=1 ligand+ macrophage
    near <- rowSums(dm[, lig_pos, drop = FALSE] <= 0.05) > 0
    rec_expr <- d[rc, ste]
    wt <- tryCatch(wilcox.test(rec_expr[near], rec_expr[!near])$p.value, error = function(e) NA_real_)
    res_all[[paste(slide, lr)]] <- data.frame(
      slide = slide, lr = lr, ligand = lg, receptor = rc,
      n_lig_pos_mac = length(lig_pos), frac_mac_lig = mean(lig_expr > 0),
      n_ste_near = sum(near), n_ste_far = sum(!near),
      rec_mean_near = mean(rec_expr[near]), rec_mean_far = mean(rec_expr[!near]),
      log2fc = log2((mean(rec_expr[near]) + 0.01) / (mean(rec_expr[!near]) + 0.01)),
      p = wt)
  }
}
res <- bind_rows(res_all)
res$fdr <- p.adjust(res$p, "BH")
print(res %>% arrange(fdr))
write.csv(res, file.path(OUT, "cosmx_lr_colocalization.csv"), row.names = FALSE)
cat("DONE\n")
