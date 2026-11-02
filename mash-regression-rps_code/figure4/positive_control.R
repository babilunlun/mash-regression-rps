#!/usr/bin/env Rscript
# positive_control.R — eQTLGen blood eQTL -> MASLD outcomes for known MASLD genes
# HSD17B13: LoF protective => higher expression => RISK, expect b>0
# MBOAT7:   lower expression => risk => higher expr => protective, expect b<0
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table)
source("/workspace/fig4/mr_core.R")
RAW <- "/workspace/fig4/raw"; QC <- "/workspace/fig4/qc"
PANEL <- file.path(RAW, "EUR38s")

pvar <- fread(file.path(RAW, "EUR38s.pvar"), skip = "#CHROM")
setnames(pvar, c("chr", "pos", "id", "ref", "alt"))
pvar[, chr := as.integer(fifelse(chr == "X", "23", chr))]
afreq <- fread(file.path(RAW, "EUR38s.afreq"), skip = "#CHROM")
setnames(afreq, c("chr", "id", "ref", "alt", "prov_ref", "alt_freq", "obs_ct"))
pvar <- merge(pvar[, .(chr, pos, id, ref, alt)], afreq[, .(id, alt_freq)], by = "id", all.x = TRUE)
setkey(pvar, chr, pos)
match_panel <- function(dt) {
  dt[, `:=`(panel_id = NA_character_, panel_eaf = NA_real_)]
  for (i in seq_len(nrow(dt))) {
    ea <- dt$ea[i]; oa <- dt$oa[i]
    if (is.na(dt$chr[i]) || is.na(dt$pos[i]) || is.na(ea) || is.na(oa)) next
    cand <- pvar[.(dt$chr[i], dt$pos[i]), nomatch = 0]
    for (j in seq_len(nrow(cand))) {
      a1 <- cand$ref[j]; a2 <- cand$alt[j]
      if ((ea == a1 && oa == a2) || (ea == a2 && oa == a1)) {
        dt$panel_id[i] <- cand$id[j]
        dt$panel_eaf[i] <- if (ea == a2) cand$alt_freq[j] else 1 - cand$alt_freq[j]
        break
      }
    }
  }
  dt
}

genes <- list(HSD17B13 = list(chr = 4L, lo = 86300000L, hi = 88500000L, expect = ">0"),
              MBOAT7   = list(chr = 19L, lo = 53000000L, hi = 56000000L, expect = "<0"))

eq <- fread(cmd = paste("zcat", shQuote(file.path(RAW, "eqtlgen_cis_fdr005.txt.gz"))))
setnames(eq, c("p", "snp", "chr", "pos", "ea", "oa", "z", "ens", "gene", "gchr", "gpos", "ncoh", "n", "fdr", "bonf"))
eq <- eq[gene %in% names(genes) & p < 5e-8]
eq[, `:=`(chr = as.integer(chr), pos = as.integer(pos), n = as.integer(n))]
bed <- eq[, .(paste0("chr", chr), pos - 1L, pos, .I)]
fwrite(bed, file.path(QC, "pc_lo.bed"), sep = "\t", col.names = FALSE)
system2(file.path(RAW, "liftOver"), c(file.path(QC, "pc_lo.bed"),
        file.path(RAW, "hg19ToHg38.over.chain.gz"),
        file.path(QC, "pc_lo38.bed"), file.path(QC, "pc_unmapped.txt")))
lo <- fread(file.path(QC, "pc_lo38.bed"), header = FALSE, col.names = c("chr38", "start", "end", "id"))
eq <- eq[lo$id]
eq[, `:=`(chr38 = as.integer(sub("^X$", "23", sub("chr", "", lo$chr38))), pos38 = lo$end)]
eq[, chr := chr38][, pos := pos38]
eq <- match_panel(eq)
eq[!is.na(panel_eaf), `:=`(se = 1 / sqrt(2 * panel_eaf * (1 - panel_eaf) * (n + z^2)),
        beta = z / sqrt(2 * panel_eaf * (1 - panel_eaf) * (n + z^2)))]
eq <- eq[!is.na(beta) & se > 0 & !is.na(panel_id)]
cat("instruments pre-clump:", nrow(eq), "\n")

clumped <- list()
for (g in names(genes)) {
  d <- eq[gene == g]
  tf <- tempfile(tmpdir = file.path(QC, "clump"))
  fwrite(unique(d[order(p), .(ID = panel_id, P = p)], by = "ID"), paste0(tf, ".txt"), sep = "\t")
  system2("plink2", c("--pfile", PANEL, "--clump", paste0(tf, ".txt"),
                      "--clump-p1", 5e-8, "--clump-p2", 1, "--clump-r2", 0.1,
                      "--clump-kb", 500, "--out", tf, "--silent"), stdout = FALSE, stderr = FALSE)
  cf <- paste0(tf, ".clumps")
  clumped[[g]] <- if (file.exists(cf)) d[panel_id %in% fread(cf)$ID] else d[which.min(p)]
  unlink(paste0(tf, "*"))
}
inst <- rbindlist(clumped)
cat("clumped:", nrow(inst), "\n")

outcomes <- list(
  finngen_fibrochir = list(f = "finngen_R12_K11_FIBROCHIRLIV.gz", n = 487866, type = "fg"),
  finngen_nafld     = list(f = "finngen_R12_NAFLD.gz", n = 487866, type = "fg"),
  gcst90091033_nafld_emerge = list(f = "GCST90091033.h.tsv.gz", n = 778614, type = "gc"),
  gcst90267352_pdff = list(f = "GCST90267352.h.tsv.gz", n = 33235, type = "gc"),
  gcst90590317_pdff = list(f = "GCST90590317.h.tsv.gz", n = 37589, type = "gc"),
  chen_pred_masld   = list(f = "pc_chen_pred_masld", n = 371405, type = "pc"),
  chen_pred_pdff    = list(f = "pc_chen_pred_pdff", n = 36675, type = "pc"))

pal <- c("A/T", "T/A", "C/G", "G/C")
res_all <- list()
for (g in names(genes)) {
  gi <- genes[[g]]
  for (oc in names(outcomes)) {
    o <- outcomes[[oc]]
    if (o$type == "fg") {
      d <- fread(cmd = sprintf("zcat %s | awk 'NR==1 || ($1==%d && $2 > %d && $2 < %d)'",
                               shQuote(file.path(RAW, o$f)), gi$chr, gi$lo, gi$hi))
      oy <- data.table(pos = d$pos, o_ea = d$alt, o_oa = d$ref, o_beta = d$beta, o_se = d$sebeta, o_p = d$pval)
    } else if (o$type == "gc") {
      d <- fread(cmd = sprintf("zcat %s | awk 'NR==1 || ($1==%d && $2 > %d && $2 < %d)'",
                               shQuote(file.path(RAW, o$f)), gi$chr, gi$lo, gi$hi))
      pick <- function(...) { nms <- c(...); nms[nms %in% names(d)][1] }
      oy <- data.table(pos = as.integer(d[[pick("hm_pos", "base_pair_location")]]),
                       o_ea = toupper(d[[pick("hm_effect_allele", "effect_allele")]]),
                       o_oa = toupper(d[[pick("hm_other_allele", "other_allele")]]),
                       o_beta = as.numeric(d[[pick("hm_beta", "beta")]]),
                       o_se = d$standard_error, o_p = as.numeric(d$p_value))
    } else {
      d <- fread(sprintf("/workspace/fig4/regions/%s.tsv", o$f))
      d <- d[gene == g]
      oy <- data.table(pos = d$pos, o_ea = d$ea, o_oa = d$oa, o_beta = d$beta, o_se = d$se, o_p = d$p)
    }
    oy[, k := paste(gi$chr, pos)]
    ix <- inst[gene == g][, k := paste(chr38, pos38)]
    m <- merge(ix, oy, by = "k")
    m <- m[!(paste(ea, oa, sep = "/") %in% pal)]
    m[, flip := (ea == o_oa & oa == o_ea)]
    m <- m[(ea == o_ea & oa == o_oa) | flip]
    if (nrow(m) == 0) { cat(sprintf("%s x %s: no overlap\n", g, oc)); next }
    m[, by := fifelse(flip, -o_beta, o_beta)]
    res <- run_mr(m$beta, m$by, m$se, m$o_se, n_x = rep(31684, nrow(m)), n_y = rep(o$n, nrow(m)),
                  eaf = m$panel_eaf)
    cat(sprintf("%-9s x %-24s n_iv=%2d  b=%+.4f se=%.4f p=%.2e  expect %s\n",
                g, oc, nrow(m), res$primary$b, res$primary$se, res$primary$p, gi$expect))
    res_all[[paste(g, oc)]] <- data.table(gene = g, outcome = oc, nsnp = nrow(m),
                                          b = res$primary$b, se = res$primary$se, p = res$primary$p,
                                          expect = gi$expect)
  }
}
fwrite(rbindlist(res_all), file.path(QC, "positive_control_results.csv"))
cat("POSITIVE_CONTROL_DONE\n")
