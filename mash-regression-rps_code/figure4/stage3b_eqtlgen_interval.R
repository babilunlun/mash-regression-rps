#!/usr/bin/env Rscript
# stage3b_eqtlgen_interval.R — eQTLGen region extraction + liftOver (FIXED join) + INTERVAL regions
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table)
RAW <- "/workspace/fig4/raw"; QC <- "/workspace/fig4/qc"; REXP <- "/workspace/fig4/regions_exp"
gmap <- fread(file.path(QC, "gene_map_gtex_liver.csv"))

## ---------- eQTLGen: probes file (gene -> hg19 chr/pos from FDR005 file) ----------
eq <- fread(cmd = paste("zcat", shQuote(file.path(RAW, "eqtlgen_cis_fdr005.txt.gz"))))
setnames(eq, c("p", "snp", "chr", "pos", "ea", "oa", "z", "ens", "gene", "gchr", "gpos", "ncoh", "n", "fdr", "bonf"))
probes <- unique(eq[gene %in% gmap$gene, .(gene, gchr, gpos)], by = "gene")
# if multiple probes per gene, keep the one with most significant min p
best <- eq[gene %in% gmap$gene, .(minp = min(p)), by = .(gene, gchr, gpos)]
probes <- best[, .SD[which.min(minp)], by = "gene"][, .(gene, gchr, gpos)]
fwrite(probes, file.path(QC, "eqtlgen_probes.tsv"), sep = "\t", col.names = FALSE)
cat("probes:", nrow(probes), "\n")

## ---------- awk extraction (GRCh37) ----------
awk_script <- '
BEGIN { FS = "\t"; OFS = "\t"
  while ((getline line < ARGV[1]) > 0) { split(line, a, "\t"); probe[a[1]] = a[2] SUBSEP a[3] }
  delete ARGV[1]
  print "gene", "chr", "pos", "rsid", "ea", "oa", "z", "p", "n" }
NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
{ g = $h["GeneSymbol"]
  if (!(g in probe)) next
  split(probe[g], pp, SUBSEP)
  c = $h["SNPChr"] + 0
  if (c != pp[1]) next
  pos = $h["SNPPos"] + 0
  d = pos - pp[2]; if (d < 0) d = -d
  if (d > 500000) next
  print g, c, pos, $h["SNP"], $h["AssessedAllele"], $h["OtherAllele"], $h["Zscore"], $h["Pvalue"], $h["NrSamples"] }'
writeLines(awk_script, file.path(RAW, "extract_eqtlgen_regions.awk"))
system(sprintf("zcat %s | awk -f %s %s > %s",
               shQuote(file.path(RAW, "eqtlgen_cis_full.txt.gz")),
               shQuote(file.path(RAW, "extract_eqtlgen_regions.awk")),
               shQuote(file.path(QC, "eqtlgen_probes.tsv")),
               shQuote(file.path(REXP, "eQTLGen_blood_eQTL__regions.tsv"))))
d <- fread(file.path(REXP, "eQTLGen_blood_eQTL__regions.tsv"))
cat("eQTLGen regions (hg19):", nrow(d), "\n")

## ---------- liftOver with CORRECT join (id = data row index) ----------
d[, rowid := .I]
fwrite(d[, .(paste0("chr", chr), pos - 1L, pos, rowid)], "/tmp/eq_reg_hg19.bed",
       sep = "\t", col.names = FALSE)
system2(file.path(RAW, "liftOver"), c("/tmp/eq_reg_hg19.bed",
        file.path(RAW, "hg19ToHg38.over.chain.gz"), "/tmp/eq_reg_hg38.bed", "/tmp/eq_reg_unmap.bed"))
lift <- fread("/tmp/eq_reg_hg38.bed", header = FALSE, col.names = c("chrs", "s", "pos38", "rowid"))
lift[, chr38 := as.integer(sub("^X$", "23", sub("chr", "", chrs)))]
d <- merge(d, lift[, .(rowid, chr38, pos38)], by = "rowid")
d <- d[!is.na(chr38)]
# sanity: CSTB must be chr21 only
cat("CSTB chrs:", paste(unique(d[gene == "CSTB"]$chr38), collapse = ","), "\n")

## ---------- panel EAF match ----------
pvar <- fread(file.path(RAW, "EUR38s.pvar"), skip = "#CHROM")
setnames(pvar, c("chr", "pos", "id", "ref", "alt"))
pvar[, chr := as.integer(fifelse(chr == "X", "23", chr))]
afreq <- fread(file.path(RAW, "EUR38s.afreq"), skip = "#CHROM")
setnames(afreq, c("chr", "id", "ref", "alt", "prov_ref", "alt_freq", "obs_ct"))
pvar <- merge(pvar[, .(chr, pos, id, ref, alt)], afreq[, .(id, alt_freq)], by = "id", all.x = TRUE)
setkey(pvar, chr, pos)
setkey(d, chr38, pos38)
m <- pvar[d, nomatch = 0]
m[, panel_eaf := fcase(ea == alt & oa == ref, alt_freq,
                       ea == ref & oa == alt, 1 - alt_freq, default = NA_real_)]
m <- m[!is.na(panel_eaf)]
out <- m[, .(gene, chr38 = chr, pos38 = pos, rsid, ea, oa, z, p, n, panel_eaf)]
fwrite(out, file.path(REXP, "eQTLGen_blood_eQTL__regions_hg38.tsv"), sep = "\t")
cat("eQTLGen regions hg38 + EAF:", nrow(out), "\n")

## ---------- INTERVAL per-gene ±500kb region files ----------
imap <- fread(file.path(RAW, "interval_core_gene_map.csv"))
for (g in unique(imap$gene)) {
  gi <- gmap[gene == g]
  if (nrow(gi) == 0) next
  accs <- strsplit(imap[gene == g, accessions], ";")[[1]]
  for (acc in accs) {
    f <- file.path(RAW, "interval", paste0(acc, ".h.tsv.gz"))
    if (!file.exists(f)) next
    d <- try(fread(cmd = paste("zcat", shQuote(f))), silent = TRUE)
    if (inherits(d, "try-error") || nrow(d) == 0) next
    d <- d[as.integer(hm_chrom) == gi$chr23 & hm_pos > gi$tss - 5e5 & hm_pos < gi$tss + 5e5]
    if (nrow(d) == 0) next
    out <- d[, .(chr = as.integer(hm_chrom), pos = hm_pos, rsid = hm_rsid,
                 ea = toupper(hm_effect_allele), oa = toupper(hm_other_allele),
                 beta = hm_beta, se = standard_error, p = as.numeric(p_value),
                 eaf = hm_effect_allele_frequency)]
    out <- out[!is.na(beta) & !is.na(se) & se > 0]
    fwrite(out, file.path(REXP, paste0("INTERVAL_pQTL:", acc, "__", g, ".tsv")), sep = "\t")
  }
}
cat("INTERVAL region files:", length(list.files(REXP, pattern = "^INTERVAL")), "\n")
cat("STAGE3B_DONE\n")
