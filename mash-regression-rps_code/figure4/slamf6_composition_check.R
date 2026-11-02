#!/usr/bin/env Rscript
# slamf6_composition_check.R — is SLAMF6 bulk decline during regression composition-driven?
# Mouse GSE305330: per-cell Slamf6 in lymphocytes, peak vs 2wk vs 6wk regression (2 models)
# Human GSE298719: per-cell SLAMF6 in T/NK, healthy vs MASLD (6 vs 7 donors)
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages({library(data.table); library(Matrix)})
OUT <- "/workspace/fig4/qc"

read_sample <- function(dir, px) {
  feat <- read.delim(gzfile(file.path(dir, paste0(px, "_features.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)
  mat <- readMM(gzfile(file.path(dir, paste0(px, "_matrix.mtx.gz"))))
  colnames(mat) <- paste0(px, "_cell", seq_len(ncol(mat)))
  rownames(mat) <- if (ncol(feat) >= 2) make.unique(feat[[2]]) else make.unique(feat[[1]])
  mat
}

# lineage markers (mouse symbols; human = uppercase)
mk <- list(
  T_cell = c("Cd3d", "Cd3e", "Cd3g"),
  NK = c("Nkg7", "Gzma", "Klrb1c"),
  B_cell = c("Cd79a", "Cd79b", "Ms4a1"),
  Myeloid = c("Lyz2", "C1qa", "Csfr1"),
  Hepatocyte = c("Alb", "Ttr", "Serpina1c"))

analyze <- function(mat, sample, model, tp, gene = "Slamf6") {
  lib <- colSums(mat)
  keep <- lib >= 500
  mat <- mat[, keep, drop = FALSE]; lib <- lib[keep]
  norm <- t(t(mat) / lib * 1e4)
  logn <- log1p(norm)
  gn <- rownames(mat)
  # lineage scores = mean log-normalized marker expression
  scores <- sapply(mk, function(g) {
    g2 <- intersect(g, gn)
    if (!length(g2)) return(rep(NA_real_, ncol(mat)))
    colMeans(logn[g2, , drop = FALSE])
  })
  scores <- as.data.table(scores)
  scores[, lineage := c("T_cell", "NK", "B_cell", "Myeloid", "Hepatocyte")[max.col(cbind(T_cell, NK, B_cell, Myeloid, Hepatocyte))]]
  lym <- scores$T_cell > 0.3 | scores$NK > 0.3 | scores$B_cell > 0.3
  slam <- as.numeric(logn[gene, ])
  slam_raw <- as.numeric(mat[gene, ])
  is_lym_lineage <- scores$lineage %in% c("T_cell", "NK", "B_cell")
  data.table(sample = sample, model = model, tp = tp,
             n_cells = ncol(mat),
             frac_lymphocyte_gate = mean(lym),
             frac_lym_lineage = mean(is_lym_lineage),
             slamf6_mean_lym = mean(slam[is_lym_lineage]),
             slamf6_pct_lym = mean(slam_raw[is_lym_lineage] > 0) * 100,
             slamf6_mean_myeloid = mean(slam[scores$lineage == "Myeloid"]),
             slamf6_mean_all = mean(slam),
             wilcox_p_lym = NA_real_)
}

mdir <- "/workspace/fig3_data/GSE305330"
mmap <- data.table(
  px = c("GSM9168811_KL70-CCl6wk", "GSM9168812_KL70-CCl8wk", "GSM9168813_KL70-CCl12wk",
         "GSM9168814_KL70-FAT8wk", "GSM9168815_KL70-FAT12wk", "GSM9168816_KL70-FAT6wk"),
  model = c("CCl4", "CCl4", "CCl4", "FAT-MASH", "FAT-MASH", "FAT-MASH"),
  tp = c("peak", "reg2wk", "reg6wk", "peak", "reg2wk", "reg6wk"))

res <- list(); cells_all <- list()
for (i in seq_len(nrow(mmap))) {
  mat <- read_sample(mdir, mmap$px[i])
  res[[i]] <- analyze(mat, mmap$px[i], mmap$model[i], mmap$tp[i])
  # keep per-cell values for tests/plots
  lib <- colSums(mat); keep <- lib >= 500
  mat <- mat[, keep, drop = FALSE]; lib <- lib[keep]
  logn <- log1p(t(t(mat) / lib * 1e4))
  gn <- rownames(mat)
  sc <- sapply(mk, function(g) colMeans(logn[intersect(g, gn), , drop = FALSE]))
  lin <- c("T_cell", "NK", "B_cell", "Myeloid", "Hepatocyte")[max.col(sc)]
  cells_all[[i]] <- data.table(sample = mmap$px[i], model = mmap$model[i], tp = mmap$tp[i],
                               lineage = lin, slamf6 = as.numeric(logn["Slamf6", ]),
                               slamf6_pos = as.numeric(mat["Slamf6", ]) > 0)
}
mres <- rbindlist(res)
cells <- rbindlist(cells_all)

# Wilcoxon peak vs reg6wk within model, lymphoid lineages only
for (mo in c("CCl4", "FAT-MASH")) {
  a <- cells[model == mo & tp == "peak" & lineage %in% c("T_cell", "NK", "B_cell")]$slamf6
  b <- cells[model == mo & tp == "reg6wk" & lineage %in% c("T_cell", "NK", "B_cell")]$slamf6
  mres[model == mo & tp == "reg6wk", wilcox_p_lym := wilcox.test(a, b)$p.value]
}
cat("=== MOUSE per-sample summary ===\n")
print(mres[, .(model, tp, n_cells, frac_lym_lineage = round(frac_lym_lineage, 3),
               slamf6_mean_lym = round(slamf6_mean_lym, 4), slamf6_pct_lym = round(slamf6_pct_lym, 2),
               slamf6_mean_all = round(slamf6_mean_all, 4), wilcox_p_lym = signif(wilcox_p_lym, 2))])

# per-lineage breakdown
cat("\n=== MOUSE Slamf6 mean by lineage x timepoint ===\n")
print(cells[, .(mean = round(mean(slamf6), 4), pct_pos = round(mean(slamf6_pos) * 100, 1), n = .N),
            by = .(model, tp, lineage)][order(model, lineage, tp)])

fwrite(mres, file.path(OUT, "slamf6_mouse_composition.csv"))
fwrite(cells, file.path(OUT, "slamf6_mouse_percell.csv"))
cat("MOUSE_DONE\n")
