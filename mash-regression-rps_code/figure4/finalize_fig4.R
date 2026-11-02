#!/usr/bin/env Rscript
# finalize_fig4.R — concordance + power + sc-eQTL annotation + PheWAS input (p-fix baked in)
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table)
QC <- "/workspace/fig4/qc"

mr <- fread(file.path(QC, "mr_results.csv"))
comp <- fread("/mnt/results/Figure2/tables/components.csv")
core <- fread("/mnt/results/Figure2/tables/core_genes_strict.txt", header = FALSE)$V1
comp_core <- comp[gene %in% core & !is.na(regr_logFC)]

m <- merge(mr, comp_core[, .(gene, regr_logFC, regr_fdr, component)], by = "gene", all.x = TRUE)
m[, `:=`(mr_sign = sign(b), regr_sign = sign(regr_logFC))]
m[, concordant := !is.na(regr_sign) & mr_sign == -regr_sign]
m[, modality := fcase(concordant == TRUE & regr_logFC < 0, "inhibitor",
                      concordant == TRUE & regr_logFC > 0, "agonist",
                      default = "discordant_or_NA")]
m[, mde80 := 2.802 * se]
m[, powered := abs(b) >= mde80]
fwrite(m, file.path(QC, "mr_results_annotated.csv"))
cat("=== MR annotated:", nrow(m), "tests |", sum(m$p < 0.05, na.rm = TRUE), "nominal |",
    sum(m$fdr < 0.05, na.rm = TRUE), "FDR<0.05 ===\n")
print(m[p < 0.05, .(n = .N, concordant = sum(concordant == TRUE, na.rm = TRUE)),
        by = .(layer_group, outcome)])

## ---- sc-eQTL cell-type annotation (GSE289173) ----
sc <- fread(cmd = "zcat /workspace/fig4/raw/gse289173/significant_sc_eQTLs_with_annotations.txt.gz")
sc_genes <- unique(sc$gene)
hits <- m[p < 0.05]
fin_genes <- unique(hits$gene)
ann <- data.table(gene = fin_genes)
for (ct in c("hepatocyte", "cholangiocyte", "stellate_cell", "endothelial_cell")) {
  col <- paste0("is_liver_eQTL_", ct)
  icol <- paste0("is_ieQTL_", ct)
  g_eq <- unique(sc[get(col) == TRUE, gene])
  g_ie <- unique(sc[get(icol) == TRUE, gene])
  ann[, paste0("sceqtl_", ct) := gene %in% g_eq]
  ann[, paste0("ieqtl_", ct) := gene %in% g_ie]
}
ann[, any_sceqtl := gene %in% sc_genes]
fwrite(ann, file.path(QC, "sceqtl_annotation.csv"))
cat("=== sc-eQTL annotation for", length(fin_genes), "nominal-hit genes ===\n")
print(ann[, .(genes_with_sceQTL = sum(any_sceqtl),
              hep = sum(sceqtl_hepatocyte), chol = sum(sceqtl_cholangiocyte),
              stel = sum(sceqtl_stellate_cell), endo = sum(sceqtl_endothelial_cell))])

## ---- PheWAS input ----
instc <- fread(file.path(QC, "instruments_clumped.csv"))
harm <- fread(file.path(QC, "harmonised_snps.csv"))
mrp <- m[p < 0.05, .(gene, layer, outcome, mr_p = p)]
lead <- merge(harm, mrp, by = c("gene", "layer", "outcome"))
lead <- lead[order(px)][, .SD[1], by = .(gene, layer, outcome)]
lead <- merge(lead, instc[, .(panel_id, rsid)], by.x = "snp", by.y = "panel_id", all.x = TRUE)
fwrite(unique(lead[, .(gene, layer, outcome, snp, rsid, chr, pos)]),
       file.path(QC, "phewas_input.tsv"), sep = "\t")
cat("phewas input:", nrow(unique(lead[, .(gene, layer, outcome)])), "pairs,",
    sum(!is.na(lead$rsid)), "with rsid\n")
cat("FINALIZE_DONE\n")
