#!/usr/bin/env Rscript
## decontx_human.R v2 — A4: DecontX ambient-RNA correction of human sc (GSE298719)
## Memory-safe redesign: everything computed PER SAMPLE inside the loop (UCell is
## per-cell rank-based => per-sample scoring identical to merged). No merged matrix.
suppressMessages({library(celda); library(Matrix); library(data.table); library(UCell)})
OUT <- "/workspace/decontx"; dir.create(OUT, showWarnings = FALSE)
DATA <- "/workspace/gse298719/GSE298719"

meta <- read.csv("/workspace/human_sc_meta.csv", stringsAsFactors = FALSE)
cat("meta:", nrow(meta), "cells |", length(unique(meta$sample)), "samples\n")
emg <- trimws(readLines("/mnt/shared-workspace/shared/human_rebuild/genes_emergent_human.txt")); emg <- emg[nzchar(emg)]
rev <- trimws(readLines("/mnt/shared-workspace/shared/human_rebuild/genes_reversal_human.txt")); rev <- rev[nzchar(rev)]

samples <- unique(meta$sample)
f_cells <- file.path(OUT, "percell_scores.csv")      # incremental append
f_cont  <- file.path(OUT, "decontx_contamination_metrics.csv")
f_stat  <- file.path(OUT, "decontx_stats_per_sample.csv")
f_det   <- file.path(OUT, "tp53inp1_detection_by_celltype_corrected.csv")
for (f in c(f_cells, f_cont, f_stat, f_det)) if (file.exists(f)) file.remove(f)

for (s in samples) {
  cat("\n=== ", s, " ===\n")
  fx <- file.path(DATA, paste0(s, "_features.tsv.gz")); bx <- file.path(DATA, paste0(s, "_barcodes.tsv.gz"))
  mx <- file.path(DATA, paste0(s, "_matrix.mtx.gz"))
  if (!file.exists(mx)) { cat("MISSING, skip\n"); next }
  feats <- read.delim(gzfile(fx), header = FALSE, stringsAsFactors = FALSE)
  bc <- read.delim(gzfile(bx), header = FALSE, stringsAsFactors = FALSE)[[1]]
  mat <- readMM(gzfile(mx)); rownames(mat) <- make.unique(feats[[2]]); colnames(mat) <- bc
  cm <- meta[meta$sample == s, ]
  keep <- intersect(colnames(mat), cm$barcode_raw)
  mat <- mat[, keep]
  cl <- cm$celltype_label[match(keep, cm$barcode_raw)]
  grp <- cm$group[match(keep, cm$barcode_raw)]
  cat("cells:", ncol(mat), "\n")
  res <- tryCatch(decontX(x = mat, z = cl), error = function(e) { cat("decontX failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  out <- if (is.list(res) && "decontXcounts" %in% names(res)) res$decontXcounts else decontXcounts(res)
  contamination <- res$contamination
  out <- round(out)
  if (!is(out, "sparseMatrix")) out <- Matrix(out, sparse = TRUE)
  cat("median contamination:", round(median(contamination), 4), "| mean:", round(mean(contamination), 4), "\n")
  rm(res); gc(verbose = FALSE)

  ## (i) contamination metrics pre/post
  nonhep_idx <- which(!grepl("Hepatocyte|Cholangiocyte", cl))
  hsc_idx <- which(cl == "HSC/myofibroblast")
  frac <- function(mm, idx, g) if (length(idx) > 0 && g %in% rownames(mm)) mean(mm[g, idx] > 0) else NA_real_
  row_cont <- data.table(
    sample = s,
    ALB_nonhep_pre = frac(mat, nonhep_idx, "ALB"),  ALB_nonhep_post = frac(out, nonhep_idx, "ALB"),
    ALB_HSC_pre = frac(mat, hsc_idx, "ALB"),        ALB_HSC_post = frac(out, hsc_idx, "ALB"),
    APOC3_HSC_pre = frac(mat, hsc_idx, "APOC3"),    APOC3_HSC_post = frac(out, hsc_idx, "APOC3"),
    TTR_HSC_pre = frac(mat, hsc_idx, "TTR"),        TTR_HSC_post = frac(out, hsc_idx, "TTR"))
  fwrite(row_cont, f_cont, append = file.exists(f_cont))
  row_stat <- data.table(sample = s, group = unique(grp), n_cells = ncol(mat),
                         median_contamination = median(contamination), mean_contamination = mean(contamination))
  fwrite(row_stat, f_stat, append = file.exists(f_stat))

  ## (ii) TP53INP1 detection by cell type (corrected counts)
  det <- aggregate(list(n = rep(1, length(cl)), det = as.numeric(out["TP53INP1", ] > 0)), by = list(celltype = cl), sum)
  det$sample <- s
  fwrite(det, f_det, append = file.exists(f_det))

  ## normalize + UCell scores (per-cell => per-sample equivalent to merged)
  libsz <- colSums(out); lnorm <- log1p(t(t(out) / libsz * 1e4))
  feats_l <- list(emergent = intersect(emg, rownames(lnorm)),
                  reversal = intersect(rev, rownames(lnorm)),
                  emergentLO_TP53INP1 = setdiff(intersect(emg, rownames(lnorm)), "TP53INP1"))
  sc <- ScoreSignatures_UCell(lnorm, features = feats_l)
  cells <- data.table(sample = s, barcode = keep, celltype = cl, group = grp,
                      emergent_c = as.numeric(sc[, 1]), reversal_c = as.numeric(sc[, 2]),
                      emergentLO_c = as.numeric(sc[, 3]),
                      TP53INP1_expr = as.numeric(lnorm["TP53INP1", ]))
  fwrite(cells, f_cells, append = file.exists(f_cells))
  rm(out, lnorm, sc, cells, mat); gc(verbose = FALSE)
}

## ---------- combine per-cell results ----------
allc <- fread(f_cells)
cat("\ncombined per-cell:", nrow(allc), "\n")

## (ii) detection table
det_all <- fread(f_det)[, .(n = sum(n), det = sum(det)), by = celltype]
det_all[, TP53INP1_pct := round(100 * det / n, 2)]
print(det_all[order(-n)])

## (iii) SAMac program MASLD vs Healthy (corrected)
samac <- allc[grepl("SAMac", celltype)]
cat("\nSAMac cells:", nrow(samac), "\n")
w_e <- wilcox.test(emergent_c ~ group, samac); w_r <- wilcox.test(reversal_c ~ group, samac)
cat(sprintf("SAMac emergent MASLD vs Healthy (corrected): median %.4f vs %.4f, P = %.3g\n",
            median(samac$emergent_c[samac$group == "MASLD"]), median(samac$emergent_c[samac$group == "Healthy"]), w_e$p.value))
cat(sprintf("SAMac reversal MASLD vs Healthy (corrected): median %.4f vs %.4f, P = %.3g\n",
            median(samac$reversal_c[samac$group == "MASLD"]), median(samac$reversal_c[samac$group == "Healthy"]), w_r$p.value))

## (iv) TP53INP1 vs emergent-LOO coupling (corrected)
cor_one <- function(x, y) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 20 || sd(x) == 0 || sd(y) == 0) return(c(rho = NA, p = NA, n = length(x)))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)); c(rho = unname(ct$estimate), p = ct$p.value, n = length(x)) }
myeloid_re <- "Kupffer|Monocyte|Dendritic|Neutrophil|SAMac|Macro|cDC|pDC"
strata <- list(all = rep(TRUE, nrow(allc)),
               myeloid = grepl(myeloid_re, allc$celltype),
               SAMac = grepl("SAMac", allc$celltype),
               HSC = allc$celltype == "HSC/myofibroblast")
coup <- do.call(rbind, lapply(names(strata), function(nm) {
  idx <- strata[[nm]] & !is.na(allc$celltype)
  r <- cor_one(allc$TP53INP1_expr[idx], allc$emergentLO_c[idx])
  data.table(stratum = nm, rho = round(r["rho"], 3), p = signif(r["p"], 3), n = r["n"])
}))
cat("\n== corrected TP53INP1 vs emergent-LOO coupling ==\n"); print(coup)

fwrite(det_all, file.path(OUT, "tp53inp1_detection_summary_corrected.csv"))
fwrite(coup, file.path(OUT, "tp53inp1_emergentLO_coupling_corrected.csv"))
fwrite(samac[, .(sample, group, emergent_c, reversal_c)], file.path(OUT, "samac_scores_corrected.csv"))
cat("DECONTX_DONE\n")
