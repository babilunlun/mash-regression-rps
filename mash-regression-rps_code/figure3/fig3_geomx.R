#!/usr/bin/env Rscript
# Figure 3 panel e (human part): GeoMx WTA ROI-level component scoring, MASH vs Healthy
suppressMessages({ library(readxl); library(dplyr); library(ggplot2) })
OUT <- "/workspace/fig3_out"
xl <- "/workspace/fig3_data/GSE287826_All_Data_WTA.xlsx"

## metadata from series matrix
sm <- readLines(gzfile("/workspace/fig3_data/GSE287826_series_matrix.txt.gz"))
get_row <- function(pat) gsub('"', "", strsplit(sm[grepl(pat, sm)][1], "\t")[[1]][-1])
gsm <- get_row("^!Sample_geo_accession")
lib <- get_row("^!Sample_description")
geno <- get_row("genotype:")
lib <- sub("Library name: ", "", lib)
meta <- data.frame(gsm, lib, genotype = geno, stringsAsFactors = FALSE) %>%
  filter(grepl("CK8", lib)) %>%
  mutate(roi = as.integer(sub(".*\\|([0-9]+)\\|.*", "\\1", lib)))
cat("mapped ROIs:", nrow(meta), " table:", paste(names(table(meta$genotype)), table(meta$genotype)), "\n")

## expression (Q3 norm)
q3 <- read_excel(xl, sheet = "Q3 Norm")
genes <- q3$TargetName
expr <- as.matrix(q3[, -1]); rownames(expr) <- genes
colnames(expr) <- gsub(" ", "", colnames(expr))  # "TMA_NASH|001|CK8_18+"
# map columns to genotype
col_roi <- as.integer(sub(".*\\|([0-9]+)\\|.*", "\\1", colnames(expr)))
meta$col <- paste0("TMA_NASH|", sprintf("%03d", meta$roi), "|CK8_18+")
stopifnot(all(meta$col %in% colnames(expr)))
expr <- expr[, meta$col]
meta$group <- meta$genotype
cat("expr:", nrow(expr), "x", ncol(expr), "\n")

## component scores: mean z-score of signature genes
rev_g <- intersect(readLines(file.path(OUT, "genes_reversal_human.txt")), rownames(expr))
eme_g <- intersect(readLines(file.path(OUT, "genes_emergent_human.txt")), rownames(expr))
cat("genes found: rev", length(rev_g), "eme", length(eme_g), "\n")
lexpr <- log2(expr + 1)
z <- t(scale(t(lexpr)))  # z per gene across ROIs
score <- data.frame(
  roi = meta$roi, group = meta$group,
  reversal = colMeans(z[rev_g, , drop = FALSE], na.rm = TRUE),
  emergent = colMeans(z[eme_g, , drop = FALSE], na.rm = TRUE))
write.csv(score, file.path(OUT, "geomx_roi_scores.csv"), row.names = FALSE)

wt_r <- wilcox.test(reversal ~ group, score)
wt_e <- wilcox.test(emergent ~ group, score)
cat("reversal MASH vs Healthy: medians", median(score$reversal[score$group == "MASH"]),
    median(score$reversal[score$group == "Healthy"]), " p =", wt_r$p.value, "\n")
cat("emergent MASH vs Healthy: medians", median(score$emergent[score$group == "MASH"]),
    median(score$emergent[score$group == "Healthy"]), " p =", wt_e$p.value, "\n")

## also: SAMac marker score (OLR1/TREM2/CD9/SPP1) as cellular-context readout
samac <- intersect(c("OLR1", "TREM2", "CD9", "SPP1", "GPNMB", "LGALS3"), rownames(expr))
score$samac <- colMeans(z[samac, , drop = FALSE], na.rm = TRUE)
hsc <- intersect(c("COL1A1", "COL3A1", "ACTA2"), rownames(expr))
score$hsc <- colMeans(z[hsc, , drop = FALSE], na.rm = TRUE)
write.csv(score, file.path(OUT, "geomx_roi_scores.csv"), row.names = FALSE)
cat("cor(reversal, samac):", cor(score$reversal, score$samac, method = "spearman"), "\n")
cat("cor(reversal, hsc):", cor(score$reversal, score$hsc, method = "spearman"), "\n")
cat("DONE\n")
