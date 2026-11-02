#!/usr/bin/env Rscript
# make_fig4_v2.R — Figure 4 v2: forest faceted by outcome scale, GTEx power-gap annotation,
# sc-eQTL layer supplement (S5), meta-analysis/evidence-tier supplement (S6)
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table); library(ggplot2); library(patchwork)
source("/workspace/fig1_theme.R")
QC <- "/workspace/fig4/qc"
OUT <- "/workspace/fig4/figs"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

layer_cols <- c(GTEx_liver_eQTL = "#DCA7EB", eQTLGen_blood_eQTL = "#A7C0DF",
                INTERVAL_pQTL = "#BFDA8E", sc_eQTL = "#E2EA4D")
oc_lab <- c(finngen_fibrochir = "FinnGen fibrosis/cirrhosis", finngen_nafld = "FinnGen MASLD",
            gcst90091033_nafld_emerge = "Ghodsian MASLD", gcst90267352_pdff = "PDFF (GCST90267352)",
            gcst90590317_pdff = "PDFF (GCST90590317)", chen_pred_masld = "Chen pred-MASLD",
            chen_pred_pdff = "Chen pred-PDFF")
oc_type <- c(finngen_fibrochir = "Case-control (log OR)", finngen_nafld = "Case-control (log OR)",
             gcst90091033_nafld_emerge = "Case-control (log OR)",
             gcst90267352_pdff = "Quantitative (beta)", gcst90590317_pdff = "Quantitative (beta)",
             chen_pred_masld = "Quantitative (beta)", chen_pred_pdff = "Quantitative (beta)")

m <- fread(file.path(QC, "mr_results_annotated.csv"))
coloc <- fread(file.path(QC, "coloc_results.csv"))
pw <- fread(file.path(QC, "phewas_associations.csv"))
meta <- fread(file.path(QC, "evidence_tiers.csv"))
scm <- fread(file.path(QC, "sceqtl_mr_results.csv"))
sus <- if (file.exists(file.path(QC, "susie_coloc_h3_recheck.csv"))) fread(file.path(QC, "susie_coloc_h3_recheck.csv")) else data.table()

## ---------- panel a: funnel + GTEx power-gap annotation ----------
instc <- fread(file.path(QC, "instruments_clumped.csv"))
inst <- fread(file.path(QC, "instruments_unclumped.csv"))
n_gtex <- nrow(instc[layer_group == "GTEx_liver_eQTL"])
funnel <- data.table(
  stage = c("Core regression genes", "Mapped to annotation", "Genes with cis instrument",
            "Clumped instruments", "Gene x outcome tests", "Nominal p<0.05", "FDR<0.05"),
  n = c(162, nrow(fread(file.path(QC, "gene_map_gtex_liver.csv"))),
        length(unique(inst$gene)), nrow(instc), nrow(m),
        sum(m$p < 0.05, na.rm = TRUE), sum(m$fdr < 0.05, na.rm = TRUE)))
funnel[, stage := factor(stage, levels = rev(stage))]
pa <- ggplot(funnel, aes(x = n, y = stage)) +
  geom_col(fill = "#A7C0DF", alpha = 0.75, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.2, family = "Liberation Sans") +
  annotate("text", x = max(funnel$n) * 0.45, y = 2.2,
           label = sprintf("GTEx liver: %d/%s instruments (%.1f%%)", n_gtex,
                           format(nrow(instc), big.mark = ","), 100 * n_gtex / nrow(instc)),
           size = 3, family = "Liberation Sans", color = "grey30", hjust = 0) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Count", y = NULL) + theme_ipf_lf()

## ---------- panel b: MR forest, faceted by outcome scale ----------
fin <- m[fdr < 0.1 | (p < 0.01 & concordant == TRUE)]
if (nrow(fin) > 30) fin <- fin[order(p)][1:30]
fin[, `:=`(lab = paste0(gene, " | ", layer_group, " | ", oc_lab[outcome]),
           otype = oc_type[outcome])]
fin <- fin[order(otype, b)]
fin[, lab := factor(lab, levels = unique(lab))]
pb <- ggplot(fin, aes(x = b, y = lab, color = layer_group)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.3, linewidth = 0.5) +
  geom_point(aes(shape = concordant == TRUE), size = 2, alpha = 0.75) +
  facet_grid(otype ~ ., scales = "free", space = "free_y") +
  scale_color_manual(values = layer_cols, name = "Evidence layer") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Direction concordant") +
  labs(x = "MR effect (outcome scale)", y = NULL) + theme_ipf_lf() +
  theme(axis.text.y = element_text(size = 5.5), strip.text = element_text(size = 7))

## ---------- panel c: coloc, H3 loci annotated with SuSiE verdict ----------
cc <- merge(m[, .(gene, layer, outcome, p, b)], coloc, by = c("gene", "layer", "outcome"))
cc[, lg := sub(":.*", "", layer)]
if (nrow(sus)) {
  cc <- merge(cc, sus[, .(gene, layer, outcome, susie_verdict = verdict)],
              by = c("gene", "layer", "outcome"), all.x = TRUE)
  cc[, pt_shape := fcase(PPH4 > 0.8, "H4",
                         !is.na(susie_verdict) & susie_verdict == "colocalised_signal_pair", "H3_rescued",
                         PPH3 > 0.8, "H3_confirmed",
                         default = "other")]
} else {
  cc[, pt_shape := fcase(PPH4 > 0.8, "H4", PPH3 > 0.8, "H3", default = "other")]
}
shp <- c(H4 = 16, H3_rescued = 17, H3_confirmed = 1, H3 = 1, other = 4)
pc <- ggplot(cc, aes(x = -log10(p), y = PPH4, color = lg)) +
  geom_hline(yintercept = 0.8, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_point(aes(shape = pt_shape), size = 2, alpha = 0.75) +
  ggrepel::geom_text_repel(data = cc[PPH4 > 0.5], aes(label = gene), size = 2.8,
                           family = "Liberation Sans", show.legend = FALSE,
                           max.overlaps = 20, box.padding = 0.4, point.padding = 0.2,
                           min.segment.length = 0.2, seed = 7) +
  scale_color_manual(values = layer_cols, name = "Evidence layer") +
  scale_shape_manual(values = shp, name = "coloc class",
                     labels = c(H4 = "PPH4>0.8", H3_rescued = "H3, SuSiE colocalised",
                                H3_confirmed = "H3, SuSiE distinct", other = "other")) +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(x = "MR -log10(p)", y = "coloc PP.H4") + theme_ipf_lf()

## ---------- panel d: PheWAS heatmap (lead-SNP level noted) ----------
pw[, trait_cat := fcase(
  grepl("liver|hepat|steatos|cirrho|fibro|MASLD|MASH|NAFLD", trait, ignore.case = TRUE), "Liver",
  grepl("leukocyte|lymphocyte|neutrophil|monocyte|eosinophil|basophil|platelet|erythrocyte|blood|immune|inflam", trait, ignore.case = TRUE), "Blood/immune",
  grepl("lipid|cholesterol|triglycer|lipoprotein|apolipo", trait, ignore.case = TRUE), "Lipids",
  grepl("diabet|glucose|insulin|glyc|HbA1c|body mass|BMI|weight|obesity|waist", trait, ignore.case = TRUE), "Metabolic",
  default = "Other")]
pws <- pw[, .(minp = min(as.numeric(p), na.rm = TRUE)), by = .(gene, trait_cat)]
pws[, ml := pmin(-log10(minp), 50)]
pd <- ggplot(pws, aes(x = trait_cat, y = gene, fill = ml)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#FAF9F3", high = "#DCA7EB", name = "-log10(min p)") +
  labs(x = "PheWAS trait category (lead-SNP level)", y = NULL) + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 6.5))

## ---------- panel e: power / MDE ----------
pe <- ggplot(m[p < 0.1], aes(x = abs(mde80), y = abs(b), color = layer_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_color_manual(values = layer_cols, name = "Evidence layer") +
  labs(x = "MDE (80% power)", y = "|MR beta|") + theme_ipf_lf()

fig4 <- (pa | pb) / (pc | pd | pe) + plot_annotation(tag_levels = "a")
save_panel(fig4, file.path(OUT, "fig4_main_v2"), width = 14, height = 10)

## ---------- S5: sc-eQTL layer ----------
# (a) FN3K blood-vs-liver eQTL direction for shared SNPs
sc <- fread(file.path(QC, "sceqtl_instruments_clumped.csv"))
eq <- fread("/workspace/fig4/regions_exp/eQTLGen_blood_eQTL__regions_hg38.tsv")
fn3k_sc <- sc[gene == "FN3K"]
eqf <- eq[gene == "FN3K"]
cmp <- list()
for (i in seq_len(nrow(fn3k_sc))) {
  ix <- fn3k_sc[i]
  h <- eqf[chr38 == ix$chr & pos38 == ix$pos]
  if (!nrow(h)) next
  h <- h[1]
  same <- (h$ea == ix$ea & h$oa == ix$oa); flip <- (h$ea == ix$oa & h$oa == ix$ea)
  if (!same && !flip) next
  cmp[[length(cmp) + 1]] <- data.table(ct = sub("sc_eQTL_", "", ix$layer), pos = ix$pos,
                                       liver_beta = ix$beta, blood_z = if (same) h$z else -h$z)
}
cmp <- rbindlist(cmp)
s5a <- ggplot(cmp, aes(x = blood_z, y = liver_beta, color = ct)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey70", linewidth = 0.4) +
  geom_point(size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c(hepatocyte = "#DCA7EB", endothelial_cell = "#E2EA4D"), name = "Liver cell type") +
  labs(x = "Blood eQTL z (eQTLGen, aligned allele)", y = "Liver sc-eQTL beta (same allele)") +
  theme_ipf_lf()

# (b) sc-eQTL MR forest (top hits)
scm[, `:=`(lab = paste0(gene, " | ", sub("sc_eQTL_", "", layer), " | ", oc_lab[outcome]),
           otype = oc_type[outcome])]
scm <- scm[order(p)]
scmt <- scm[1:min(12, .N)]
scmt <- scmt[order(b)]
scmt[, lab := factor(lab, levels = unique(lab))]
s5b <- ggplot(scmt, aes(x = b, y = lab)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.3, linewidth = 0.5, color = "#E2EA4D") +
  geom_point(size = 2, alpha = 0.75, color = "#E2EA4D") +
  labs(x = "sc-eQTL MR effect (outcome scale)", y = NULL) + theme_ipf_lf() +
  theme(axis.text.y = element_text(size = 6))
s5 <- (s5a | s5b) + plot_annotation(tag_levels = "a")
save_panel(s5, file.path(OUT, "fig4S5_sceqtl_layer"), width = 11, height = 4.5)

## ---------- S6: meta-analysis / evidence tiers ----------
mt <- meta[!is.na(clin_p) & !is.na(disc_p)]
mt[, tier_short := sub("_.*", "", tier)]
tier_cols <- c(T1 = "#E2EA4D", T2 = "#BFDA8E", T3 = "#A7C0DF", T4 = "#DCA7EB", T5 = "#C8CACB", T6 = "#C7B0C7")
s6a <- ggplot(mt, aes(x = -log10(disc_p), y = -log10(clin_p), color = tier_short)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey70", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = 3, color = "grey70", linewidth = 0.4) +
  geom_point(size = 2, alpha = 0.75) +
  ggrepel::geom_text_repel(data = mt[tier_short %in% c("T1", "T2") | clin_p < 1e-3],
                           aes(label = gene), size = 2.8, family = "Liberation Sans",
                           show.legend = FALSE, max.overlaps = 20, seed = 7) +
  scale_color_manual(values = tier_cols, name = "Tier") +
  labs(x = "Discovery -log10(p) (Chen pred-MASLD)", y = "Clinical MASLD meta -log10(p)") +
  theme_ipf_lf()
td <- meta[, .N, by = tier][order(tier)]
td[, tier := factor(tier, levels = rev(tier))]
s6b <- ggplot(td, aes(x = N, y = tier)) +
  geom_col(fill = "#BEF3C5", alpha = 0.75, width = 0.7) +
  geom_text(aes(label = N), hjust = -0.15, size = 3.2, family = "Liberation Sans") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Gene x layer combinations", y = NULL) + theme_ipf_lf() +
  theme(axis.text.y = element_text(size = 6.5))
s6 <- (s6a | s6b) + plot_annotation(tag_levels = "a")
save_panel(s6, file.path(OUT, "fig4S6_meta_tiers"), width = 11, height = 4.5)

cat("FIG4_V2_DONE\n")
