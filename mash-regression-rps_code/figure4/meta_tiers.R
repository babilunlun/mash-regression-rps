#!/usr/bin/env Rscript
# meta_tiers.R — cross-outcome meta-analysis (Stouffer, sqrt-n weighted) + evidence tiers
# Clinical MASLD set: finngen_fibrochir, finngen_nafld, gcst90091033 (case-control cohorts)
# Discovery: chen_pred_masld (ML-predicted). PDFF set: gcst90267352, gcst90590317, chen_pred_pdff.
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table)
QC <- "/workspace/fig4/qc"

m <- fread(file.path(QC, "mr_results_annotated.csv"))
cc <- fread(file.path(QC, "coloc_results.csv"))
sc <- if (file.exists(file.path(QC, "sceqtl_mr_results.csv"))) fread(file.path(QC, "sceqtl_mr_results.csv")) else data.table()

N <- c(finngen_fibrochir = 487866, finngen_nafld = 487866, gcst90091033_nafld_emerge = 778614,
       gcst90267352_pdff = 33235, gcst90590317_pdff = 37589,
       chen_pred_masld = 371405, chen_pred_pdff = 36675)
clinical <- c("finngen_fibrochir", "finngen_nafld", "gcst90091033_nafld_emerge")
discovery <- "chen_pred_masld"
pdff <- c("gcst90267352_pdff", "gcst90590317_pdff", "chen_pred_pdff")

stouffer <- function(b, se, n) {
  z <- b / se
  w <- sqrt(n)
  zc <- sum(w * z) / sqrt(sum(w^2))
  2 * pnorm(-abs(zc))
}

meta <- m[, {
  cl <- .SD[outcome %in% clinical]
  pd <- .SD[outcome %in% pdff]
  di <- .SD[outcome %in% discovery]
  list(
    clin_p = if (nrow(cl) >= 2) stouffer(cl$b, cl$se, N[cl$outcome]) else NA_real_,
    clin_dir = if (nrow(cl)) mean(sign(cl$b) == sign(cl$b[which.min(cl$p)])) else NA_real_,
    clin_minp = if (nrow(cl)) min(cl$p) else NA_real_,
    pdff_p = if (nrow(pd) >= 2) stouffer(pd$b, pd$se, N[pd$outcome]) else NA_real_,
    pdff_dir = if (nrow(pd)) mean(sign(pd$b) == sign(pd$b[which.min(pd$p)])) else NA_real_,
    disc_p = if (nrow(di)) di$p[1] else NA_real_,
    disc_b = if (nrow(di)) di$b[1] else NA_real_,
    best_fdr = min(fdr, na.rm = TRUE),
    best_outcome = outcome[which.min(p)],
    best_b = b[which.min(p)],
    concordant_any = any(concordant == TRUE, na.rm = TRUE),
    modality = modality[which.min(p)]
  )
}, by = .(gene, layer_group)]

ccb <- cc[, .(PPH4_max = max(PPH4), PPH3_max = max(PPH3),
              best_coloc_outcome = outcome[which.max(PPH4)]), by = .(gene, layer_group = sub(":.*", "", layer))]
meta <- merge(meta, ccb, by = c("gene", "layer_group"), all.x = TRUE)

# SuSiE-coloc recheck: rescue colocalised signal pairs, confirm distinct-signal loci
susf <- file.path(QC, "susie_coloc_h3_recheck.csv")
if (file.exists(susf) && file.size(susf) > 10) {
  sus <- fread(susf)
  sus[, layer_group := sub(":.*", "", layer)]
  susv <- sus[, .(susie_verdict = verdict[which.max(ifelse(is.na(max_PPH4), -1, max_PPH4))],
                  susie_PPH4 = max(max_PPH4, na.rm = TRUE)), by = .(gene, layer_group)]
  meta <- merge(meta, susv, by = c("gene", "layer_group"), all.x = TRUE)
  meta[!is.na(susie_verdict) & susie_verdict == "colocalised_signal_pair", PPH4_max := susie_PPH4]
}

if (nrow(sc)) {
  scs <- sc[, .(sceqtl_nominal = any(p < 0.05),
                sceqtl_best_p = min(p),
                sceqtl_ct = paste(unique(sub("sc_eQTL_", "", layer)[p < 0.05]), collapse = ";")),
            by = gene]
  meta <- merge(meta, scs, by = "gene", all.x = TRUE)
}

meta[, tier := fcase(
  best_fdr < 0.05 & concordant_any & PPH4_max > 0.8, "T1_FDR+concordant+coloc",
  best_fdr < 0.05 & PPH4_max > 0.8, "T2_FDR+coloc_dirNA",
  best_fdr < 0.05 & concordant_any & PPH4_max > 0.5, "T2_FDR+concordant+coloc05",
  best_fdr < 0.05 & PPH3_max > 0.8, "T4_FDR_but_LDconfounded",
  best_fdr < 0.05 & concordant_any, "T3_FDR+concordant",
  best_fdr < 0.05, "T5_FDR_only",
  default = "T6_nominal_or_null")]
meta <- meta[order(tier, best_fdr)]
fwrite(meta, file.path(QC, "evidence_tiers.csv"))

cat("=== tier distribution ===\n")
print(meta[, .N, by = tier][order(tier)])
cat("\n=== T1/T2 targets ===\n")
cols <- c("gene", "layer_group", "best_outcome", "best_b", "best_fdr", "clin_p", "clin_dir",
          "PPH4_max", "modality", "tier", "sceqtl_ct")
cols <- cols[cols %in% names(meta)]
print(meta[tier %like% "T1|T2", ..cols])
cat("\n=== key genes clinical-MASLD meta (discovery excluded) ===\n")
print(meta[gene %in% c("FN3K", "TP53INP1", "SLAMF6", "ICAM1", "CTSS", "TSPAN3", "HAVCR2", "CD52", "OAS2"),
           .(gene, layer_group, clin_p = signif(clin_p, 2), clin_dir = round(clin_dir, 2),
             disc_p = signif(disc_p, 2), best_fdr = signif(best_fdr, 2), tier)][order(clin_p)])
cat("META_TIERS_DONE\n")
