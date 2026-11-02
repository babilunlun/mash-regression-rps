#!/usr/bin/env Rscript
# stage5_proxy_mr.R — LD proxy rescue + harmonisation + MR (7 outcomes)
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table)
source("/workspace/fig4/mr_core.R")
RAW <- "/workspace/fig4/raw"; QC <- "/workspace/fig4/qc"; REG <- "/workspace/fig4/regions"
REXP <- "/workspace/fig4/regions_exp"; PANEL <- file.path(RAW, "EUR38s")
dir.create(file.path(QC, "ld"), showWarnings = FALSE)

instc <- fread(file.path(QC, "instruments_clumped.csv"))
if (!"layer_group" %in% names(instc)) instc[, layer_group := sub(":.*", "", layer)]
gmap <- fread(file.path(QC, "gene_map_gtex_liver.csv"))
gmap[, chr23 := as.integer(sub("^X$", "23", sub("chr", "", chr)))]

## ---- exposure pools ----
sp_all <- fread(file.path(QC, "gtex_liver_core_pairs_all.csv"))
pv <- strsplit(sp_all$variant_id, "_", fixed = TRUE)
gtex_pool <- data.table(gene = gmap$gene[match(sub("\\..*", "", sp_all$gene_id), sub("\\..*", "", gmap$ens))],
                        chr = as.integer(sub("^X$", "23", sub("chr", "", sapply(pv, `[`, 1)))),
                        pos = as.integer(sapply(pv, `[`, 2)),
                        oa = sapply(pv, `[`, 3), ea = sapply(pv, `[`, 4),
                        beta = sp_all$slope, se = sp_all$slope_se, p = sp_all$pval_nominal,
                        eaf = sp_all$maf, layer_group = "GTEx_liver_eQTL")
eq_pool <- fread(file.path(REXP, "eQTLGen_blood_eQTL__regions_hg38.tsv"))
eq_pool[, `:=`(se = 1 / sqrt(2 * panel_eaf * (1 - panel_eaf) * (n + z^2)),
               beta = z / sqrt(2 * panel_eaf * (1 - panel_eaf) * (n + z^2)))]
eq_pool <- eq_pool[!is.na(beta) & se > 0,
                   .(gene, chr = chr38, pos = pos38, ea, oa, beta, se, p, eaf = panel_eaf,
                     layer_group = "eQTLGen_blood_eQTL")]
pq_pool <- list()
for (f in list.files(REXP, pattern = "^INTERVAL.*\\.tsv$")) {
  d <- fread(file.path(REXP, f))
  acc <- sub("__.*", "", sub("INTERVAL_pQTL:", "", f))
  gene <- sub(".*__", "", sub("\\.tsv$", "", f))
  d <- d[!is.na(beta) & !is.na(se) & se > 0]
  pq_pool[[f]] <- d[, .(gene, chr = as.integer(chr), pos, ea, oa, beta, se, p, eaf,
                        layer_group = paste0("INTERVAL_pQTL:", acc))]
}
pq_pool <- rbindlist(pq_pool)
exposure_pool <- rbindlist(list(gtex_pool, eq_pool, pq_pool), fill = TRUE)
exposure_pool[, k := paste(chr, pos)]
cat("exposure pool:", nrow(exposure_pool), "rows\n")

## ---- outcome variant sets ----
outcomes <- c("finngen_fibrochir", "finngen_nafld", "gcst90091033_nafld_emerge",
              "gcst90267352_pdff", "gcst90590317_pdff", "chen_pred_masld", "chen_pred_pdff")
oy_sets <- list()
for (oc in outcomes) {
  d <- fread(file.path(REG, paste0(oc, ".tsv")))
  d[, k := paste(chr, pos)]
  oy_sets[[oc]] <- d
}

## ---- LD cache: instruments vs panel variants in +-1Mb ----
ld_cache <- list()
for (g in unique(instc$gene)) {
  gi <- gmap[gene == g][1]
  inst_g <- instc[gene == g]
  pool_g <- exposure_pool[gene == g]
  if (nrow(pool_g) == 0 || is.na(gi$chr23)) next
  tf <- tempfile(tmpdir = file.path(QC, "ld"))
  fwrite(data.table(ids = unique(inst_g$panel_id)), paste0(tf, "_inst.txt"), col.names = FALSE)
  system2("plink2", c("--pfile", PANEL,
                      "--chr", gi$chr23, "--from-bp", max(0, gi$tss - 1e6), "--to-bp", gi$tss + 1e6,
                      "--r2-unphased", "--ld-snp-list", paste0(tf, "_inst.txt"),
                      "--ld-window", 1000000, "--ld-window-r2", 0.5,
                      "--out", tf, "--silent"), stdout = FALSE, stderr = FALSE)
  vf <- paste0(tf, ".vcor")
  if (file.exists(vf)) {
    v <- fread(vf)
    ld_cache[[g]] <- v[, .(inst_id = ID_A, chr_b = CHROM_B, pos_b = POS_B, r2 = UNPHASED_R2)]
  }
  unlink(paste0(tf, "*"))
}
ld_all <- rbindlist(ld_cache)
cat("LD pairs:", nrow(ld_all), "\n")
fwrite(ld_all, file.path(QC, "instrument_ld_cache.csv"))

## ---- harmonise with proxy rescue ----
pal <- c("A/T", "T/A", "C/G", "G/C")
results <- list()
for (oc in outcomes) {
  oy_all <- oy_sets[[oc]]
  oy_by_gene <- split(oy_all, oy_all$gene)
  for (i in seq_len(nrow(instc))) {
    ix <- instc[i]
    if (is.na(ix$chr) || is.na(ix$pos)) next
    og <- oy_by_gene[[ix$gene]]
    if (is.null(og)) next
    pool_gl <- exposure_pool[gene == ix$gene & layer_group == ix$layer]
    if (nrow(pool_gl) == 0) next
    ldg <- ld_all[inst_id == ix$panel_id]
    cand <- data.table()
    if (any(!is.na(pool_gl$chr) & pool_gl$chr == ix$chr & pool_gl$pos == ix$pos))
      cand <- rbind(cand, pool_gl[chr == ix$chr & pos == ix$pos][, `:=`(r2 = 1)])
    if (nrow(ldg)) {
      pr <- merge(pool_gl, ldg[r2 >= 0.8, .(chr_b, pos_b, r2)],
                  by.x = c("chr", "pos"), by.y = c("chr_b", "pos_b"))
      cand <- rbind(cand, pr, fill = TRUE)
    }
    if (nrow(cand) == 0) next
    cand <- cand[k %in% og$k]
    if (nrow(cand) == 0) next
    setorder(cand, -r2, p)
    for (j in seq_len(nrow(cand))) {
      cx <- cand[j]
      if (paste(cx$ea, cx$oa, sep = "/") %in% pal) next
      ocand <- og[k == paste(cx$chr, cx$pos)]
      if (nrow(ocand) == 0) next
      hit <- NA; flip <- FALSE
      for (jj in seq_len(nrow(ocand))) {
        if (ocand$ea[jj] == cx$ea && ocand$oa[jj] == cx$oa) { hit <- jj; flip <- FALSE; break }
        if (ocand$ea[jj] == cx$oa && ocand$oa[jj] == cx$ea) { hit <- jj; flip <- TRUE; break }
      }
      if (is.na(hit)) next
      by <- ocand$beta[hit]; sey <- ocand$se[hit]
      if (is.na(by) || is.na(sey) || sey <= 0) next
      if (flip) by <- -by
      results[[length(results) + 1]] <- data.table(
        gene = ix$gene, layer = ix$layer, layer_group = ix$layer_group, tier = ix$tier,
        outcome = oc, snp = ix$panel_id, proxy_k = paste(cx$chr, cx$pos), proxy_r2 = cx$r2,
        chr = ix$chr, pos = ix$pos, ea = cx$ea, oa = cx$oa,
        bx = cx$beta, bxse = cx$se, px = cx$p,
        by = by, byse = sey, py = ocand$p[hit], eaf_x = cx$eaf, eaf_y = ocand$eaf[hit])
      break
    }
  }
}
harm <- rbindlist(results, fill = TRUE)
cat("harmonised rows:", nrow(harm), "| proxy usage:", round(mean(harm$proxy_r2 < 1) * 100, 1), "%\n")
fwrite(harm, file.path(QC, "harmonised_snps.csv"))

## ---- MR ----
outcome_n <- list(finngen_fibrochir = 487866, finngen_nafld = 487866,
                  gcst90091033_nafld_emerge = 778614,
                  gcst90267352_pdff = 33235, gcst90590317_pdff = 37589,
                  chen_pred_masld = 371405, chen_pred_pdff = 36675)
mr_rows <- list()
groups <- unique(harm[, .(gene, layer, layer_group, tier, outcome)])
for (i in seq_len(nrow(groups))) {
  g <- groups[i]
  d <- harm[gene == g$gene & layer == g$layer & outcome == g$outcome]
  d <- d[!is.na(bx) & !is.na(by) & !is.na(bxse) & !is.na(byse) & bxse > 0 & byse > 0 & bx != 0]
  if (nrow(d) == 0) next
  nx <- instc[gene == g$gene & layer == g$layer]$n[1]
  res <- run_mr(d$bx, d$by, d$bxse, d$byse,
                n_x = rep(nx, nrow(d)), n_y = rep(outcome_n[[g$outcome]], nrow(d)), eaf = d$eaf_x)
  row <- data.table(gene = g$gene, layer = g$layer, layer_group = g$layer_group,
                    tier = g$tier, outcome = g$outcome, nsnp = res$nsnp,
                    mean_proxy_r2 = mean(d$proxy_r2),
                    b = res$primary$b, se = res$primary$se, p = res$primary$p,
                    method = res$primary$method)
  if (!is.null(res$primary$Q)) row[, `:=`(Q = res$primary$Q, Q_df = res$primary$Q_df, Q_p = res$primary$Q_p)]
  if (!is.null(res$egger)) row[, `:=`(egger_b = res$egger$b, egger_se = res$egger$se, egger_p = res$egger$p,
                                      egger_int = res$egger$intercept, egger_int_p = res$egger$intercept_p)]
  if (!is.null(res$wmedian)) row[, `:=`(wm_b = res$wmedian$b, wm_se = res$wmedian$se, wm_p = res$wmedian$p)]
  if (!is.null(res$steiger)) row[, `:=`(steiger_rx = res$steiger$r_exposure,
                                        steiger_ry = res$steiger$r_outcome,
                                        steiger_ok = res$steiger$correct_direction)]
  mr_rows[[length(mr_rows) + 1]] <- row
}
mr <- rbindlist(mr_rows, fill = TRUE)
mr[, fdr := p.adjust(p, method = "BH"), by = .(layer_group, outcome)]
mr[, fdr_global := p.adjust(p, method = "BH")]
fwrite(mr, file.path(QC, "mr_results.csv"))
cat("=== MR results ===\n")
print(mr[, .(n_tests = .N, sig_fdr05 = sum(fdr < 0.05, na.rm = TRUE),
             sig_p05 = sum(p < 0.05, na.rm = TRUE)), by = .(layer_group, outcome)])
cat("STAGE5_DONE\n")
