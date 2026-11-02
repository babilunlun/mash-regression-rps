#!/usr/bin/env Rscript
# sc_expression_8genes.R — cell-type expression of 8 scorecard genes in human + mouse liver sc
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages({library(data.table); library(Matrix)})
OUT <- "/workspace/fig5"

GENES_H <- c("TP53INP1", "SLAMF6", "FN3K", "MYO18A", "OAS2", "TSPAN3", "CTSS", "ICAM1")
GENES_M <- c("Trp53inp1", "Slamf6", "Fn3k", "Myo18a", "Oas2", "Tspan3", "Ctss", "Icam1")

mk_h <- list(T_cell = c("CD3D", "CD3E", "CD3G"), NK = c("NKG7", "GZMA", "KLRB1"),
             B_cell = c("CD79A", "CD79B", "MS4A1"), Kupffer = c("MARCO", "VSIG4", "CD163"),
             SAMac = c("TREM2", "CD9", "OLR1", "GPNMB"), Monocyte_MF = c("LYZ", "S100A8", "S100A9", "VCAN"),
             HSC_MF = c("COL1A1", "ACTA2", "LRAT", "COL3A1"), Endothelial = c("PECAM1", "FLT1", "VWF"),
             Cholangiocyte = c("KRT19", "EPCAM", "SOX9"), Hepatocyte = c("ALB", "TTR", "SERPINA1"))
mk_m <- lapply(mk_h, function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x)))))
names(mk_m) <- names(mk_h)

read_sample <- function(dir, px, has_bc = FALSE) {
  feat <- read.delim(gzfile(file.path(dir, paste0(px, "_features.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)
  mat <- readMM(gzfile(file.path(dir, paste0(px, "_matrix.mtx.gz"))))
  colnames(mat) <- paste0(px, "_c", seq_len(ncol(mat)))
  rownames(mat) <- if (ncol(feat) >= 2) make.unique(feat[[2]]) else make.unique(feat[[1]])
  mat
}

process <- function(dir, px, group, genes, mks) {
  mat <- read_sample(dir, px)
  lib <- colSums(mat); keep <- lib >= 500
  mat <- mat[, keep, drop = FALSE]; lib <- lib[keep]
  logn <- log1p(t(t(mat) / lib * 1e4))
  gn <- rownames(mat)
  sc <- sapply(mks, function(g) { g2 <- intersect(g, gn); if (!length(g2)) return(rep(0, ncol(mat))); colMeans(logn[g2, , drop = FALSE]) })
  lin <- colnames(sc)[max.col(sc)]
  # require minimum lineage score to call
  maxsc <- apply(sc, 1, max)
  lin[maxsc < 0.15] <- "Other"
  g2 <- intersect(genes, gn)
  res <- lapply(g2, function(g) {
    v <- as.numeric(logn[g, ])
    data.table(gene = g, lineage = lin, expr = v)[, .(mean_expr = mean(expr), pct_expr = mean(expr > 0) * 100, n = .N),
                                                   by = .(gene, lineage)]
  })
  r <- rbindlist(res)
  r[, `:=`(sample = px, group = group)]
  r
}

## human
hdir <- "/workspace/fig3_data/GSE298719"
hpx <- sub("_features.tsv.gz$", "", list.files(hdir, pattern = "_features.tsv.gz$"))
hpx <- hpx[grepl("FRPBiopsy", hpx)]
hgrp <- sub("^.*_(Healthy|MASLD)$", "\\1", hpx)
hres <- rbindlist(lapply(seq_along(hpx), function(i) process(hdir, hpx[i], hgrp[i], GENES_H, mk_h)))
hres[, species := "human"]
fwrite(hres, file.path(OUT, "sc_expr_human_8genes.csv"))
cat("human done:", nrow(hres), "\n")

## mouse
mdir <- "/workspace/fig3_data/GSE305330"
mmap <- data.table(
  px = c("GSM9168811_KL70-CCl6wk", "GSM9168812_KL70-CCl8wk", "GSM9168813_KL70-CCl12wk",
         "GSM9168814_KL70-FAT8wk", "GSM9168815_KL70-FAT12wk", "GSM9168816_KL70-FAT6wk"),
  group = c("CCl4_peak", "CCl4_reg2wk", "CCl4_reg6wk", "FATMASH_peak", "FATMASH_reg2wk", "FATMASH_reg6wk"))
mres <- rbindlist(lapply(seq_len(nrow(mmap)), function(i) process(mdir, mmap$px[i], mmap$group[i], GENES_M, mk_m)))
mres[, species := "mouse"]
fwrite(mres, file.path(OUT, "sc_expr_mouse_8genes.csv"))
cat("mouse done:", nrow(mres), "\n")

## quick look at the two leads
print(hres[gene %in% c("TP53INP1", "SLAMF6") & n > 50,
           .(group, gene, lineage, mean_expr = round(mean_expr, 3), pct = round(pct_expr, 1))][order(gene, group, -mean_expr)][1:24])
cat("SC_EXPR_DONE\n")
