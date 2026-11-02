#!/usr/bin/env Rscript
# slamf6_human_check.R — SLAMF6 per-cell expression in human liver T/NK, healthy vs MASLD
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages({library(data.table); library(Matrix)})
OUT <- "/workspace/fig4/qc"
DATA <- "/workspace/fig3_data/GSE298719"

mk <- list(T_cell = c("CD3D", "CD3E", "CD3G"), NK = c("NKG7", "GZMA", "KLRB1"),
           B_cell = c("CD79A", "CD79B", "MS4A1"), Myeloid = c("LYZ", "C1QA", "CSF1R"),
           Hepatocyte = c("ALB", "TTR", "SERPINA1"))

feats <- list.files(DATA, pattern = "_features.tsv.gz$")
px <- sub("_features.tsv.gz$", "", feats)
px <- px[grepl("FRPBiopsy", px)]
meta <- data.table(sample = px,
                   group = sub("^.*_(Healthy|MASLD)$", "\\1", px))
print(meta[, .N, by = group])

per_sample <- list(); per_cell <- list()
for (s in px) {
  feat <- read.delim(gzfile(file.path(DATA, paste0(s, "_features.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)
  mat <- readMM(gzfile(file.path(DATA, paste0(s, "_matrix.mtx.gz"))))
  rownames(mat) <- if (ncol(feat) >= 2) make.unique(feat[[2]]) else make.unique(feat[[1]])
  lib <- colSums(mat); keep <- lib >= 500
  mat <- mat[, keep, drop = FALSE]; lib <- lib[keep]
  logn <- log1p(t(t(mat) / lib * 1e4))
  gn <- rownames(mat)
  sc <- sapply(mk, function(g) { g2 <- intersect(g, gn); colMeans(logn[g2, , drop = FALSE]) })
  lin <- c("T_cell", "NK", "B_cell", "Myeloid", "Hepatocyte")[max.col(sc)]
  slam <- as.numeric(logn["SLAMF6", ]); slam_raw <- as.numeric(mat["SLAMF6", ])
  grp <- meta[sample == s]$group
  is_tnk <- lin %in% c("T_cell", "NK")
  per_sample[[s]] <- data.table(sample = s, group = grp, n_cells = ncol(mat),
                                frac_TNK = mean(is_tnk),
                                slamf6_mean_TNK = mean(slam[is_tnk]),
                                slamf6_pct_TNK = mean(slam_raw[is_tnk] > 0) * 100,
                                slamf6_mean_all = mean(slam))
  per_cell[[s]] <- data.table(sample = s, group = grp, lineage = lin, slamf6 = slam)
}
ps <- rbindlist(per_sample)
pc <- rbindlist(per_cell)
fwrite(ps, file.path(OUT, "slamf6_human_per_sample.csv"))
fwrite(pc, file.path(OUT, "slamf6_human_per_cell.csv"))

cat("=== HUMAN per-donor summary ===\n")
print(ps[, .(sample, group, n_cells, frac_TNK = round(frac_TNK, 3),
             slamf6_TNK = round(slamf6_mean_TNK, 4), pct_TNK = round(slamf6_pct_TNK, 1),
             slamf6_all = round(slamf6_mean_all, 4))][order(group, sample)])
cat("\n=== donor-level tests (n=6 healthy vs 7 MASLD) ===\n")
w1 <- wilcox.test(ps[group == "Healthy"]$slamf6_mean_TNK, ps[group == "MASLD"]$slamf6_mean_TNK)
w2 <- wilcox.test(ps[group == "Healthy"]$frac_TNK, ps[group == "MASLD"]$frac_TNK)
w3 <- wilcox.test(ps[group == "Healthy"]$slamf6_mean_all, ps[group == "MASLD"]$slamf6_mean_all)
cat("per-cell SLAMF6 in T/NK:  healthy", round(median(ps[group=="Healthy"]$slamf6_mean_TNK),4),
    "vs MASLD", round(median(ps[group=="MASLD"]$slamf6_mean_TNK),4), " Wilcoxon p =", signif(w1$p.value, 3), "\n")
cat("T/NK fraction:            healthy", round(median(ps[group=="Healthy"]$frac_TNK),3),
    "vs MASLD", round(median(ps[group=="MASLD"]$frac_TNK),3), " Wilcoxon p =", signif(w2$p.value, 3), "\n")
cat("bulk-equivalent SLAMF6:   healthy", round(median(ps[group=="Healthy"]$slamf6_mean_all),4),
    "vs MASLD", round(median(ps[group=="MASLD"]$slamf6_mean_all),4), " Wilcoxon p =", signif(w3$p.value, 3), "\n")
cat("HUMAN_DONE\n")
