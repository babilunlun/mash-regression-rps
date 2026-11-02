.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages(library(data.table))
QC <- "/workspace/fig4/qc"; REXP <- "/workspace/fig4/regions_exp"

# sc-eQTL FN3K instruments
sc <- fread(file.path(QC, "sceqtl_instruments_clumped.csv"))
fn3k_sc <- sc[gene == "FN3K"]
cat("=== FN3K sc-eQTL instruments ===\n")
print(fn3k_sc[, .(layer, chr, pos, ea, oa, beta = signif(beta, 3), p = signif(p, 1))])

# eQTLGen blood: effect of the SAME SNPs on FN3K
eq <- fread(file.path(REXP, "eQTLGen_blood_eQTL__regions_hg38.tsv"))
eqf <- eq[gene == "FN3K"]
cat("\n=== same SNPs in eQTLGen blood ===\n")
for (i in seq_len(nrow(fn3k_sc))) {
  ix <- fn3k_sc[i]
  hit <- eqf[chr38 == ix$chr & pos38 == ix$pos]
  if (nrow(hit)) {
    h <- hit[1]
    # align: sc ea/oa vs eq ea/oa
    same <- (h$ea == ix$ea & h$oa == ix$oa)
    flip <- (h$ea == ix$oa & h$oa == ix$ea)
    z_al <- if (same) h$z else if (flip) -h$z else NA
    cat(ix$layer, ix$pos, "| sc beta(sc-ea) =", signif(ix$beta, 3),
        "| blood z(aligned to sc-ea) =", ifelse(is.na(z_al), "allele-mismatch", signif(z_al, 3)),
        "| blood p =", signif(h$p, 1), "\n")
  } else cat(ix$layer, ix$pos, "| not in blood region file\n")
}

# eQTLGen FN3K top instruments for reference
instc <- fread(file.path(QC, "instruments_clumped.csv"))
cat("\n=== FN3K blood instruments (top 5) ===\n")
print(instc[gene == "FN3K"][order(p)][1:5, .(chr, pos, ea, oa, p = signif(p, 1))])

# CYP1A2 in blood MR
m <- fread(file.path(QC, "mr_results_annotated.csv"))
cat("\n=== CYP1A2 blood eQTL MR ===\n")
print(m[gene == "CYP1A2", .(layer_group, outcome, nsnp, b = signif(b, 3), p = signif(p, 2))][order(p)][1:6])
