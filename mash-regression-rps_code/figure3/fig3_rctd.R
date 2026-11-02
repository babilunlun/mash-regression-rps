#!/usr/bin/env Rscript
# Figure 3 panel e (mouse part): RCTD deconvolution of Visium peak vs regression slides
# using the annotated mouse sc reference; + per-spot component scoring
suppressMessages({ library(Seurat); library(dplyr); library(spacexr); library(Matrix); library(UCell) })
OUT <- "/workspace/fig3_out"; DATA <- "/workspace/fig3_data/GSE306327"

## ---- reference ----
mo <- readRDS(file.path(OUT, "mouse_sc_scored.rds"))
mo <- JoinLayers(mo)
counts <- GetAssayData(mo, layer = "counts")
ct <- factor(mo$celltype_label)
names(ct) <- colnames(mo)
# drop rare/artifact types for reference stability
keep_t <- names(which(table(ct) >= 25))
keep_cells <- names(ct)[ct %in% keep_t]
ref <- Reference(counts[, keep_cells], ct[keep_cells], colSums(counts[, keep_cells]))
cat("reference:", length(keep_cells), "cells,", length(keep_t), "types\n")

## ---- spatial ----
read_visium <- function(px, posdir) {
  mat <- readMM(gzfile(file.path(DATA, paste0(px, "_matrix.mtx.gz"))))
  feat <- read.delim(gzfile(file.path(DATA, paste0(px, "_features.tsv.gz"))), header = FALSE)
  bc <- read.delim(gzfile(file.path(DATA, paste0(px, "_barcodes.tsv.gz"))), header = FALSE)[[1]]
  rownames(mat) <- make.unique(feat[[2]]); colnames(mat) <- bc
  pos <- read.csv(file.path(DATA, posdir, "spatial", "tissue_positions.csv"))
  rownames(pos) <- pos$barcode
  pos <- pos[colnames(mat), ]
  coords <- data.frame(x = pos$pxl_col_in_fullres, y = pos$pxl_row_in_fullres)
  rownames(coords) <- colnames(mat)
  list(counts = mat, coords = coords, pos = pos)
}
slides <- list(peak = read_visium("GSM9197415_peak-6-w-2", "peak"),
               regression = read_visium("GSM9197416_regression-6-2-w", "regression"))

rctd_out <- list()
for (sn in names(slides)) {
  sl <- slides[[sn]]
  puck <- SpatialRNA(sl$coords, sl$counts, colSums(sl$counts))
  myRCTD <- create.RCTD(puck, ref, max_cores = 8)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
  rctd_out[[sn]] <- myRCTD
  w <- myRCTD@results$weights
  write.csv(as.matrix(w), file.path(OUT, paste0("rctd_weights_", sn, ".csv")))
  cat(sn, "spots:", nrow(w), "\n")
}
saveRDS(rctd_out, file.path(OUT, "rctd_results_tmp.rds"))
file.rename(file.path(OUT, "rctd_results_tmp.rds"), file.path(OUT, "rctd_results.rds"))

## ---- per-spot component scores (UCell on spatial counts) ----
rev_g <- intersect(readLines(file.path(OUT, "genes_reversal_mouse.txt")), unique(unlist(lapply(slides, function(s) rownames(s$counts)))))
eme_g <- intersect(readLines(file.path(OUT, "genes_emergent_mouse.txt")), unique(unlist(lapply(slides, function(s) rownames(s$counts)))))
cat("spatial genes: rev", length(rev_g), "eme", length(eme_g), "\n")
for (sn in names(slides)) {
  sl <- slides[[sn]]
  d <- sl$counts
  d <- d[rowSums(d) > 0, ]
  d <- log1p(sweep(d, 2, colSums(d), "/") * 1e4)
  sc <- ScoreSignatures_UCell(d, features = list(reversal = rev_g, emergent = eme_g), name = "_prog")
  spot_df <- data.frame(barcode = colnames(d), sc, check.names = FALSE)
  write.csv(spot_df, file.path(OUT, paste0("spot_scores_", sn, ".csv")), row.names = FALSE)
}
cat("DONE\n")
