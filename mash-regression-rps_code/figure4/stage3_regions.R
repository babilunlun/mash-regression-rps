#!/usr/bin/env Rscript
# stage3_regions.R — gene map + per-gene ±1Mb region extracts for 7 outcomes (GRCh38)
# + Chen GRCh37->GRCh38 full-file liftOver + region extraction
# + eQTLGen region extraction (awk) + liftOver with CORRECT row-id join
# + INTERVAL per-gene ±500kb region files
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table)
RAW <- "/workspace/fig4/raw"; QC <- "/workspace/fig4/qc"; REG <- "/workspace/fig4/regions"
REXP <- "/workspace/fig4/regions_exp"
dir.create(REG, showWarnings = FALSE); dir.create(REXP, showWarnings = FALSE)

## ---------- gene map: core genes -> GTEx liver annotation ----------
core <- fread("/mnt/results/Figure2/tables/core_genes_strict.txt", header = FALSE)$V1
eg <- fread(cmd = paste("zcat", shQuote(file.path(RAW, "GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz"))))
aliases <- c(FYB1 = "FYB", MACROH2A2 = "H2AFY2", TCIM = "C8orf59")
eg_map <- eg[, .(gene = gene_name, ens = gene_id, chr = gene_chr, tss = gene_start)]
# TSS by strand: for - strand the TSS is gene_end
eg_map <- eg[, .(gene = gene_name, ens = gene_id, chr = gene_chr,
                 tss = fifelse(strand == "-", gene_end, gene_start))]
hit <- eg_map[match(core, eg_map$gene)]
miss <- is.na(hit$gene)
for (a in names(aliases)) {
  idx <- which(core == a)
  if (length(idx) && is.na(hit$gene[idx])) {
    h2 <- eg_map[eg_map$gene == aliases[a]]
    if (nrow(h2)) hit[idx] <- h2[1]
  }
}
gmap <- data.table(gene = core, ens = hit$ens, chr = hit$chr, tss = hit$tss)
gmap <- gmap[!is.na(ens)]
gmap[, chr23 := as.integer(sub("^X$", "23", sub("chr", "", chr)))]
cat("gene map:", nrow(gmap), "/", length(core), "\n")
fwrite(gmap, file.path(QC, "gene_map_gtex_liver.csv"))

## ---------- Chen liftOver (full files, GRCh37 -> GRCh38) ----------
chen_lift <- function(tag) {
  f <- sprintf("%s/chen_raw/Raw summary statistics/GWAS in UK Biobank/gwas_%s.tsv", RAW, tag)
  bed <- sprintf("/tmp/chen_%s_hg19.bed", tag)
  system(sprintf("awk -F'\\t' 'NR>1 {print \"chr\"$2, $3-1, $3, NR}' OFS='\\t' %s > %s", shQuote(f), bed))
  system2(file.path(RAW, "liftOver"), c(bed, file.path(RAW, "hg19ToHg38.over.chain.gz"),
          sprintf("/tmp/chen_%s_hg38.bed", tag), sprintf("/tmp/chen_%s_unmap.bed", tag)))
  d <- fread(f)
  d[, rowid := .I]
  lift <- fread(sprintf("/tmp/chen_%s_hg38.bed", tag), header = FALSE,
                col.names = c("chrs", "s", "pos38", "bid"))
  lift[, rowid := bid - 1L]  # BED id = NR = data row + 1
  lift[, chr38 := as.integer(sub("^X$", "23", sub("chr", "", chrs)))]
  m <- merge(d, lift[, .(rowid, chr38, pos38)], by = "rowid")
  m <- m[!is.na(chr38)]
  cat(tag, "lifted:", nrow(m), "| rs738409@43928847:", nrow(m[chr38 == 22 & pos38 == 43928847]), "\n")
  m
}
for (tag in c("pred", "pdff")) {
  m <- chen_lift(tag)
  oc <- if (tag == "pred") "chen_pred_masld" else "chen_pred_pdff"
  out <- list()
  for (i in seq_len(nrow(gmap))) {
    g <- gmap[i]
    sub <- m[chr38 == g$chr23 & pos38 > g$tss - 1e6 & pos38 < g$tss + 1e6]
    if (nrow(sub) == 0) next
    out[[g$gene]] <- sub[, .(gene = g$gene, chr = chr38, pos = pos38, rsid = variant_id,
                             ea = toupper(effect_allele), oa = toupper(other_allele),
                             beta, se = standard_error, p = 10^-neg_log_10_p_value,
                             eaf = effect_allele_frequency, n)]
  }
  res <- rbindlist(out)
  fwrite(res, file.path(REG, paste0(oc, ".tsv")), sep = "\t")
  cat(oc, ":", nrow(res), "rows,", length(out), "genes\n")
  # positive-control gene windows (HSD17B13, MBOAT7)
  pcg <- list(HSD17B13 = list(chr = 4L, lo = 86300000L, hi = 88500000L),
              MBOAT7 = list(chr = 19L, lo = 53000000L, hi = 56000000L))
  out2 <- list()
  for (g in names(pcg)) {
    gi <- pcg[[g]]
    sub <- m[chr38 == gi$chr & pos38 > gi$lo & pos38 < gi$hi]
    out2[[g]] <- sub[, .(gene = g, chr = chr38, pos = pos38, rsid = variant_id,
                         ea = toupper(effect_allele), oa = toupper(other_allele),
                         beta, se = standard_error, p = 10^-neg_log_10_p_value,
                         eaf = effect_allele_frequency, n)]
  }
  fwrite(rbindlist(out2), file.path(REG, paste0("pc_", oc, ".tsv")), sep = "\t")
  rm(m); gc()
}

## ---------- other 5 outcomes: direct GRCh38 region extraction ----------
extract_region <- function(d, gene, chr23, tss, ea_col, oa_col, b_col, se_col, p_col,
                           eaf_col = NULL, n_col = NULL, rsid_col = NULL) {
  sub <- d[chr == chr23 & pos > tss - 1e6 & pos < tss + 1e6]
  if (nrow(sub) == 0) return(NULL)
  data.table(gene = gene, chr = sub$chr, pos = sub$pos,
             rsid = if (!is.null(rsid_col)) sub[[rsid_col]] else NA_character_,
             ea = toupper(sub[[ea_col]]), oa = toupper(sub[[oa_col]]),
             beta = sub[[b_col]], se = sub[[se_col]], p = sub[[p_col]],
             eaf = if (!is.null(eaf_col)) sub[[eaf_col]] else NA_real_,
             n = if (!is.null(n_col)) sub[[n_col]] else NA_real_)
}

# FinnGen (GRCh38): #chrom pos ref alt rsids pval beta sebeta af_alt
for (oc in c("K11_FIBROCHIRLIV", "NAFLD")) {
  d <- fread(cmd = paste("zcat", shQuote(file.path(RAW, paste0("finngen_R12_", oc, ".gz")))),
             select = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"))
  setnames(d, c("chr", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"))
  d[, chr := as.integer(sub("^X$", "23", sub("chr", "", chr)))]
  out <- rbindlist(lapply(seq_len(nrow(gmap)), function(i) {
    g <- gmap[i]
    extract_region(d, g$gene, g$chr23, g$tss, "alt", "ref", "beta", "sebeta", "pval", "af_alt", rsid_col = "rsids")
  }))
  nm <- if (oc == "K11_FIBROCHIRLIV") "finngen_fibrochir" else "finngen_nafld"
  fwrite(out, file.path(REG, paste0(nm, ".tsv")), sep = "\t")
  cat(nm, ":", nrow(out), "\n")
  rm(d); gc()
}
# GWAS Catalog harmonised (GRCh38)
for (acc in c("GCST90091033", "GCST90267352", "GCST90590317")) {
  f_acc <- file.path(RAW, paste0(acc, ".h.tsv.gz"))
  ok <- file.exists(f_acc) && file.size(f_acc) > 1e6 &&
        system(paste("gzip -t", shQuote(f_acc)), ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  if (!ok) { cat(acc, "SKIPPED (missing or corrupt)\n"); next }
  d <- fread(cmd = paste("zcat", shQuote(f_acc)))
  pick <- function(...) { nms <- c(...); nms[nms %in% names(d)][1] }
  d <- d[, .(chr = as.integer(get(pick("hm_chrom", "chromosome"))),
             pos = as.integer(get(pick("hm_pos", "base_pair_location"))),
             rsid = get(pick("hm_rsid", "variant_id")),
             ea = toupper(get(pick("hm_effect_allele", "effect_allele"))),
             oa = toupper(get(pick("hm_other_allele", "other_allele"))),
             beta = as.numeric(get(pick("hm_beta", "beta"))),
             se = as.numeric(get(pick("standard_error"))),
             p = as.numeric(get(pick("p_value"))),
             eaf = as.numeric(get(pick("hm_effect_allele_frequency", "effect_allele_frequency"))))]
  d <- d[!is.na(chr) & !is.na(beta) & !is.na(se) & se > 0]
  out <- rbindlist(lapply(seq_len(nrow(gmap)), function(i) {
    g <- gmap[i]
    sub <- d[chr == g$chr23 & pos > g$tss - 1e6 & pos < g$tss + 1e6]
    if (nrow(sub) == 0) return(NULL)
    data.table(gene = g$gene, chr = sub$chr, pos = sub$pos, rsid = sub$rsid,
               ea = sub$ea, oa = sub$oa, beta = sub$beta, se = sub$se, p = sub$p,
               eaf = sub$eaf, n = NA_real_)
  }))
  nm <- switch(acc, GCST90091033 = "gcst90091033_nafld_emerge",
               GCST90267352 = "gcst90267352_pdff", GCST90590317 = "gcst90590317_pdff")
  fwrite(out, file.path(REG, paste0(nm, ".tsv")), sep = "\t")
  cat(nm, ":", nrow(out), "\n")
  rm(d); gc()
}
cat("STAGE3_OUTCOMES_DONE\n")
