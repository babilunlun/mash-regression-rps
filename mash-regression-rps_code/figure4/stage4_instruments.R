#!/usr/bin/env Rscript
# stage4_instruments.R — build cis instruments (3 layers) + LD clump (ALL FIXES BAKED IN)
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table)
RAW <- "/workspace/fig4/raw"; QC <- "/workspace/fig4/qc"; REXP <- "/workspace/fig4/regions_exp"
PANEL <- file.path(RAW, "EUR38s")
dir.create(file.path(QC, "clump"), showWarnings = FALSE, recursive = TRUE)

gmap <- fread(file.path(QC, "gene_map_gtex_liver.csv"))
gmap <- gmap[!is.na(chr23)]

# panel variant index (GRCh38)
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
    if (is.na(dt$chr[i]) || is.na(dt$pos[i]) || is.na(ea) || is.na(oa) || !nzchar(ea) || !nzchar(oa)) next
    cand <- pvar[.(dt$chr[i], dt$pos[i]), nomatch = 0]
    if (nrow(cand) == 0) next
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

## ---------- Layer 1: GTEx v8 Liver eQTL ----------
sp <- fread(cmd = paste("zcat", shQuote(file.path(RAW, "GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.gz"))))
sp[, ens := sub("\\..*", "", gene_id)]
gmap[, ensu := sub("\\..*", "", ens)]
sp <- sp[ens %in% gmap$ensu]
fwrite(sp, file.path(QC, "gtex_liver_core_pairs_all.csv"))
parse_vid <- function(v) {
  p <- strsplit(v, "_", fixed = TRUE)
  list(chr = sub("^X$", "23", sub("chr", "", sapply(p, `[`, 1))), pos = as.integer(sapply(p, `[`, 2)),
       ref = sapply(p, `[`, 3), alt = sapply(p, `[`, 4))
}
pv <- parse_vid(sp$variant_id)
sp[, `:=`(chr = as.integer(pv$chr), pos = pv$pos, ref = pv$ref, alt = pv$alt)]
sp[, gene := gmap$gene[match(ens, gmap$ensu)]]
strict <- sp[pval_nominal < 5e-8]
relaxed <- sp[pval_nominal < 1e-6 & !gene %in% strict$gene]
mk_gtex <- function(d, tier) {
  data.table(gene = d$gene, chr = d$chr, pos = d$pos, ea = d$alt, oa = d$ref,
             beta = d$slope, se = d$slope_se, p = d$pval_nominal, eaf = d$maf,
             n = 208L, layer = "GTEx_liver_eQTL", tier = tier)
}
gtex <- rbind(mk_gtex(strict, "strict"), mk_gtex(relaxed, "relaxed"))
gtex <- match_panel(gtex)
cat("GTEx instruments:", nrow(gtex), "| panel-matched:", sum(!is.na(gtex$panel_id)), "\n")

## ---------- Layer 2: eQTLGen blood eQTL ----------
eq <- fread(cmd = paste("zcat", shQuote(file.path(RAW, "eqtlgen_cis_fdr005.txt.gz"))))
setnames(eq, c("p", "snp", "chr", "pos", "ea", "oa", "z", "ens", "gene", "gchr", "gpos", "ncoh", "n", "fdr", "bonf"))
eq_core <- eq[gene %in% gmap$gene & p < 5e-8]
cat("eQTLGen core-gene cis pairs (5e-8):", nrow(eq_core), "across", length(unique(eq_core$gene)), "genes\n")
eqi <- data.table()
if (nrow(eq_core) > 0) {
  bed <- data.table(chr = paste0("chr", eq_core$chr), s = eq_core$pos - 1, e = eq_core$pos,
                    id = seq_len(nrow(eq_core)))
  fwrite(bed, file.path(RAW, "eqtlgen_inst_hg19.bed"), sep = "\t", col.names = FALSE)
  system2(file.path(RAW, "liftOver"),
          c(file.path(RAW, "eqtlgen_inst_hg19.bed"), file.path(RAW, "hg19ToHg38.over.chain.gz"),
            file.path(RAW, "eqtlgen_inst_hg38.bed"), file.path(RAW, "eqtlgen_inst_unmap.bed")))
  lifted <- fread(file.path(RAW, "eqtlgen_inst_hg38.bed"), header = FALSE,
                  col.names = c("chrs", "s", "pos38", "id"))
  lifted[, chr38 := as.integer(sub("^X$", "23", sub("chr", "", chrs)))]
  eq_core[, rowid := seq_len(.N)]
  eqm <- merge(eq_core, lifted[, .(rowid = id, chr38, pos38)], by = "rowid", all.x = TRUE)
  eqi <- data.table(gene = eqm$gene, chr = eqm$chr38, pos = eqm$pos38,
                    ea = eqm$ea, oa = eqm$oa, z = eqm$z, p = eqm$p, n = as.integer(eqm$n),
                    layer = "eQTLGen_blood_eQTL", tier = "strict")
  eqi <- eqi[!is.na(chr) & !is.na(pos)]
  eqi <- match_panel(eqi)
  eqi[!is.na(panel_eaf), `:=`(
    se = 1 / sqrt(2 * panel_eaf * (1 - panel_eaf) * (n + z^2)),
    beta = z / sqrt(2 * panel_eaf * (1 - panel_eaf) * (n + z^2)))]
  eqi[, eaf := panel_eaf]
  eqi[, z := NULL]
  cat("eQTLGen instruments lifted+matched:", sum(!is.na(eqi$panel_id)), "/", nrow(eqi), "\n")
}

## ---------- Layer 3: INTERVAL pQTL ----------
imap <- fread(file.path(RAW, "interval_core_gene_map.csv"))
pq_list <- list()
for (g in unique(imap$gene)) {
  accs <- strsplit(imap[gene == g, accessions], ";")[[1]]
  for (acc in accs) {
    f <- file.path(REXP, paste0("INTERVAL_pQTL:", acc, "__", g, ".tsv"))
    if (!file.exists(f)) next
    d <- try(fread(f), silent = TRUE)
    if (inherits(d, "try-error") || nrow(d) == 0) next
    cis <- d[p < 5e-8 & !is.na(beta) & !is.na(se) & se > 0]
    if (nrow(cis) == 0) next
    pq_list[[paste(g, acc)]] <- data.table(
      gene = g, chr = as.integer(cis$chr), pos = cis$pos, ea = cis$ea, oa = cis$oa,
      beta = cis$beta, se = cis$se, p = cis$p, eaf = cis$eaf, n = 3301L,
      layer = paste0("INTERVAL_pQTL:", acc), tier = "strict")
  }
}
pqi <- rbindlist(pq_list)
if (nrow(pqi)) pqi <- match_panel(pqi)
cat("INTERVAL instruments:", nrow(pqi), "| panel-matched:", sum(!is.na(pqi$panel_id)), "\n")

## ---------- combine + F stat ----------
inst <- rbindlist(list(gtex, eqi, pqi), fill = TRUE)
inst <- inst[!is.na(panel_id) & !is.na(beta) & !is.na(se) & se > 0]
inst[, f_stat := (beta / se)^2]
inst <- inst[f_stat > 10]
cat("instruments (F>10):", nrow(inst), "\n")
print(inst[, .(n_snps = .N, n_genes = length(unique(gene)), minF = min(f_stat)),
           by = .(layer_group = sub(":.*", "", layer), tier)])
fwrite(inst, file.path(QC, "instruments_unclumped.csv"))

## ---------- LD clump per gene x layer (FIXED: ID column + tab sep + dedup) ----------
clumped_list <- list()
combos <- unique(inst[, .(gene, layer, tier)])
for (i in seq_len(nrow(combos))) {
  g <- combos$gene[i]; ly <- combos$layer[i]; tier <- combos$tier[i]
  d <- inst[gene == g & layer == ly]
  if (nrow(d) == 1) { clumped_list[[paste(g, ly)]] <- d; next }
  p1 <- if (tier == "strict") 5e-8 else 1e-6
  tf <- tempfile(tmpdir = file.path(QC, "clump"))
  fwrite(unique(d[order(p), .(ID = panel_id, P = p)], by = "ID"), paste0(tf, ".txt"), sep = "\t")
  system2("plink2", c("--pfile", PANEL, "--clump", paste0(tf, ".txt"),
                      "--clump-p1", p1, "--clump-p2", 1,
                      "--clump-r2", 0.1, "--clump-kb", 500,
                      "--out", tf, "--silent"), stdout = FALSE, stderr = FALSE)
  cf <- paste0(tf, ".clumps")
  if (file.exists(cf)) {
    cl <- fread(cf)
    clumped_list[[paste(g, ly)]] <- d[panel_id %in% cl$ID]
  } else {
    clumped_list[[paste(g, ly)]] <- d[which.min(p)]
  }
  unlink(paste0(tf, "*"))
}
instc <- rbindlist(clumped_list, fill = TRUE)
instc[, layer_group := sub(":.*", "", layer)]
cat("clumped instruments:", nrow(instc), "across", length(unique(instc$gene)), "genes x",
    length(unique(instc$layer)), "layers\n")
print(instc[, .(n_snps = .N, n_genes = length(unique(gene))), by = .(layer_group, tier)])
fwrite(instc, file.path(QC, "instruments_clumped.csv"))
cat("STAGE4_DONE\n")
