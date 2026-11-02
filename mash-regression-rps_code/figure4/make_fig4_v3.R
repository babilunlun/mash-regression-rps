#!/usr/bin/env Rscript
# make_fig4_v3.R — Figure 4 v3: adds genetic-evidence scorecard as main panel f
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table); library(ggplot2); library(patchwork)
source("/workspace/fig1_theme.R")
QC <- "/workspace/fig4/qc"
OUT <- "/workspace/fig4/figs"

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
sus <- fread(file.path(QC, "susie_coloc_h3_recheck.csv"))

## ---------- panel a: funnel ----------
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
  annotate("text", x = max(funnel$n) * 0.42, y = 2.2,
           label = sprintf("GTEx liver: %d/%s (%.1f%%)", n_gtex,
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

## ---------- panel c: coloc with SuSiE verdicts ----------
cc <- merge(m[, .(gene, layer, outcome, p, b)], coloc, by = c("gene", "layer", "outcome"))
cc[, lg := sub(":.*", "", layer)]
cc <- merge(cc, sus[, .(gene, layer, outcome, susie_verdict = verdict)],
            by = c("gene", "layer", "outcome"), all.x = TRUE)
cc[, pt_shape := fcase(PPH4 > 0.8, "H4",
                       !is.na(susie_verdict) & susie_verdict == "colocalised_signal_pair", "H3_rescued",
                       PPH3 > 0.8, "H3_confirmed",
                       default = "other")]
shp <- c(H4 = 16, H3_rescued = 17, H3_confirmed = 1, other = 4)
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

## ---------- panel d: PheWAS heatmap ----------
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

## ---------- panel f: genetic evidence scorecard ----------
# dimensions: MR FDR | coloc | clinical replication | direction | tissue consistency | PheWAS clean | multi-layer
sc_genes <- c("TP53INP1", "FN3K", "SLAMF6", "MYO18A", "OAS2", "TSPAN3", "CTSS", "ICAM1")
mt8 <- meta[gene %in% sc_genes]
# best row per gene by tier then fdr
setorder(mt8, tier, best_fdr)
best <- mt8[, .SD[1], by = gene]
# multi-layer: FDR<0.05 in >=2 layers, or coloc PPH4>0.8 in >=2 layers
ml_fdr <- m[fdr < 0.05, .(n_layer = uniqueN(layer_group)), by = gene]
ml_coloc <- coloc[PPH4 > 0.8, .(n_cl = uniqueN(sub(":.*", "", layer))), by = gene]
sus_ok <- sus[verdict == "colocalised_signal_pair", unique(gene)]
# PheWAS clinical disease flags
pwc <- fread(file.path(QC, "phewas_associations_categorized.csv"))
clin_flag <- pwc[cat == "clinical/other" & !grepl("circumference|height|weight|BMI|body", trait, ignore.case = TRUE),
                 unique(gene)]
score <- CJ(gene = sc_genes, dim = c("MR FDR<0.05", "coloc PPH4>0.8", "Clinical replication",
                                     "Direction concordant", "Tissue-consistent",
                                     "PheWAS clean", "Multi-layer"))
score[, status := {
  g <- gene; d <- dim
  br <- best[gene == g]
  fdr_ok <- nrow(br) && isTRUE(br$best_fdr < 0.05)
  coloc_ok <- nrow(br) && (isTRUE(br$PPH4_max > 0.8) || g %in% sus_ok)
  clin_ok <- nrow(br) && isTRUE(br$clin_p < 0.05)
  dir_ok <- nrow(br) && isTRUE(br$modality %in% c("agonist", "inhibitor"))
  tissue <- if (g == "FN3K") "fail" else "nd"
  phewas_ok <- !(g %in% clin_flag)
  multi <- (ml_fdr[gene == g, n_layer] >= 2) || (ml_coloc[gene == g, n_cl] >= 2)
  if (length(multi) == 0) multi <- FALSE
  n_layers_tested <- uniqueN(m[gene == g]$layer_group)
  val <- switch(d,
    "MR FDR<0.05" = fdr_ok, "coloc PPH4>0.8" = coloc_ok, "Clinical replication" = clin_ok,
    "Direction concordant" = dir_ok, "Tissue-consistent" = tissue,
    "PheWAS clean" = phewas_ok, "Multi-layer" = multi)
  # test-not-applicable rules (documented in README):
  # SLAMF6 direction: bulk regression logFC invalid for lymphocyte-lineage markers (composition-confounded, see slamf6_* tables)
  if (d == "Direction concordant" && g == "SLAMF6") val <- "nd"
  # multi-layer untestable when only one layer exists for the gene
  if (d == "Multi-layer" && n_layers_tested < 2) val <- "nd"
  if (d == "Tissue-consistent") val else if (identical(val, "nd")) "nd" else ifelse(isTRUE(val), "pass", "fail")
}, by = .(gene, dim)]
score[, status := factor(status, levels = c("pass", "fail", "nd"))]
# order: core causal quartet (FDR/coloc/clinical/direction) passes first, then total passes
core_dims <- c("MR FDR<0.05", "coloc PPH4>0.8", "Clinical replication", "Direction concordant")
ord <- score[, .(core = sum(status == "pass" & dim %in% core_dims),
                 np = sum(status == "pass"), nf = sum(status == "fail")), by = gene][
       order(-core, -np, nf)]
score[, gene := factor(gene, levels = rev(ord$gene))]
score[, dim := factor(dim, levels = c("MR FDR<0.05", "coloc PPH4>0.8", "Clinical replication",
                                      "Direction concordant", "Tissue-consistent",
                                      "PheWAS clean", "Multi-layer"))]
pf <- ggplot(score, aes(x = dim, y = gene, fill = status)) +
  geom_tile(color = "white", linewidth = 0.5, alpha = 0.75) +
  geom_text(aes(label = fcase(status == "pass", "\u2713", status == "fail", "\u2717", default = "n.d.")),
            size = 3.2, family = "Liberation Sans") +
  scale_fill_manual(values = c(pass = "#BFDA8E", fail = "#DCA7EB", nd = "#C8CACB"),
                    labels = c(pass = "Pass", fail = "Fail", nd = "No data"), name = NULL) +
  labs(x = "Genetic evidence criterion", y = NULL) + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

fig4 <- (pa | pb) / (pc | pd) / (pe | pf) + plot_annotation(tag_levels = "a")
save_panel(fig4, file.path(OUT, "fig4_main_v3"), width = 13, height = 13.5)
fwrite(dcast(score, gene ~ dim, value.var = "status"), file.path(QC, "genetic_scorecard.csv"))
cat("FIG4_V3_DONE\n")
