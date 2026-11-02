#!/usr/bin/env Rscript
# add_rsids.R — recover rsIDs for clumped instruments (for PheWAS)
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table)
QC <- "/workspace/fig4/qc"; REXP <- "/workspace/fig4/regions_exp"; RAW <- "/workspace/fig4/raw"

instc <- fread(file.path(QC, "instruments_clumped.csv"))
if ("rsid" %in% names(instc)) instc[, rsid := NULL]
instc[, rsid := NA_character_]

# eQTLGen: rsid from region pool (GRCh38 pos match)
pool <- fread(file.path(REXP, "eQTLGen_blood_eQTL__regions_hg38.tsv"))
pool[, k := paste(chr38, pos38)]
eq_map <- unique(pool[!is.na(rsid), .(k, rsid)])
setkey(eq_map, k)
idx <- instc$layer_group == "eQTLGen_blood_eQTL"
instc[idx, rsid := eq_map[paste(instc$chr[idx], instc$pos[idx]), rsid]]

# INTERVAL: rsid from per-gene region files
for (f in list.files(REXP, pattern = "^INTERVAL.*\\.tsv$")) {
  d <- fread(file.path(REXP, f))
  acc <- sub("__.*", "", sub("INTERVAL_pQTL:", "", f))
  lay <- paste0("INTERVAL_pQTL:", acc)
  d[, k := paste(as.integer(chr), pos)]
  m <- unique(d[!is.na(rsid) & rsid != "", .(k, rsid)])
  setkey(m, k)
  idx <- instc$layer == lay
  instc[idx, rsid := m[paste(instc$chr[idx], instc$pos[idx]), rsid]]
}

# GTEx: rsid from egenes top-variant annotation where positions match
eg <- fread(cmd = paste("zcat", shQuote(file.path(RAW, "GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz"))))
eg <- eg[, .(gene_id, rsid_top = rs_id_dbSNP151_GRCh38p7, variant_id)]
pv <- strsplit(eg$variant_id, "_", fixed = TRUE)
eg[, `:=`(chr = as.integer(sub("^X$", "23", sub("chr", "", sapply(pv, `[`, 1)))),
          pos = as.integer(sapply(pv, `[`, 2)))]
gmap <- fread(file.path(QC, "gene_map_gtex_liver.csv"))
eg[, gene := gmap$gene[match(sub("\\..*", "", gene_id), sub("\\..*", "", gmap$ens))]]
gt <- which(instc$layer_group == "GTEx_liver_eQTL")
for (i in gt) {
  hit <- eg[gene == instc$gene[i] & chr == instc$chr[i] & pos == instc$pos[i]]
  if (nrow(hit)) instc$rsid[i] <- hit$rsid_top[1]
}

cat("rsID recovery rate:\n")
print(instc[, .(n = .N, with_rsid = sum(!is.na(rsid)),
              pct = round(mean(!is.na(rsid)) * 100, 1)), by = layer_group])
fwrite(instc, file.path(QC, "instruments_clumped.csv"))
cat("RSID_DONE\n")
