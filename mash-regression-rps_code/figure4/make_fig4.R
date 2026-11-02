#!/usr/bin/env Rscript
# make_fig4.R — Figure 4 main (5 panels) + supplements, ipf-lf-paper-rules theme
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table); library(ggplot2); library(patchwork)
source("/workspace/fig1_theme.R")
QC <- "/workspace/fig4/qc"
OUT <- "/workspace/fig4/figs"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

layer_cols <- c(GTEx_liver_eQTL = "#DCA7EB", eQTLGen_blood_eQTL = "#A7C0DF",
                INTERVAL_pQTL = "#BFDA8E")
oc_lab <- c(finngen_fibrochir = "FinnGen fibrosis/cirrhosis", finngen_nafld = "FinnGen MASLD",
            gcst90091033_nafld_emerge = "Ghodsian MASLD", gcst90267352_pdff = "PDFF (GCST90267352)",
            gcst90590317_pdff = "PDFF (GCST90590317)", chen_pred_masld = "Chen pred-MASLD",
            chen_pred_pdff = "Chen pred-PDFF")

m <- fread(file.path(QC, "mr_results_annotated.csv"))
coloc <- if (file.exists(file.path(QC, "coloc_results.csv"))) fread(file.path(QC, "coloc_results.csv")) else data.table()
pw <- if (file.exists(file.path(QC, "phewas_associations.csv"))) fread(file.path(QC, "phewas_associations.csv")) else data.table()
ann <- if (file.exists(file.path(QC, "sceqtl_annotation.csv"))) fread(file.path(QC, "sceqtl_annotation.csv")) else data.table()

## ---------- panel a: funnel ----------
instc <- fread(file.path(QC, "instruments_clumped.csv"))
inst <- fread(file.path(QC, "instruments_unclumped.csv"))
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
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Count", y = NULL) + theme_ipf_lf()

## ---------- panel b: MR forest ----------
fin <- m[fdr < 0.1 | (p < 0.01 & concordant == TRUE)]
if (nrow(fin) > 30) fin <- fin[order(p)][1:30]
fin[, lab := paste0(gene, " | ", layer_group, " | ", oc_lab[outcome])]
fin <- fin[order(b)]
fin[, lab := factor(lab, levels = unique(lab))]
pb <- ggplot(fin, aes(x = b, y = lab, color = layer_group)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.3, linewidth = 0.5) +
  geom_point(aes(shape = concordant == TRUE), size = 2, alpha = 0.75) +
  scale_color_manual(values = layer_cols, name = "Evidence layer") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Direction concordant") +
  labs(x = "MR effect (beta, outcome scale)", y = NULL) + theme_ipf_lf() +
  theme(axis.text.y = element_text(size = 5.5))

## ---------- panel c: coloc ----------
if (nrow(coloc)) {
  cc <- merge(m[, .(gene, layer, outcome, p, b)], coloc, by = c("gene", "layer", "outcome"))
  cc[, lg := sub(":.*", "", layer)]
  pc <- ggplot(cc, aes(x = -log10(p), y = PPH4, color = lg)) +
    geom_hline(yintercept = 0.8, linetype = 2, color = "grey60", linewidth = 0.4) +
    geom_point(size = 2, alpha = 0.75) +
    ggrepel::geom_text_repel(data = cc[PPH4 > 0.5], aes(label = gene), size = 2.8,
                             family = "Liberation Sans", show.legend = FALSE,
                             max.overlaps = 20, box.padding = 0.4, point.padding = 0.2,
                             min.segment.length = 0.2, seed = 7) +
    scale_color_manual(values = layer_cols, name = "Evidence layer") +
    scale_y_continuous(limits = c(0, 1.02)) +
    labs(x = "MR -log10(p)", y = "coloc PP.H4") + theme_ipf_lf()
} else {
  pc <- ggplot() + annotate("text", x = 0, y = 0, label = "No coloc pairs") + theme_ipf_lf()
}

## ---------- panel d: PheWAS heatmap ----------
if (nrow(pw)) {
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
    labs(x = "PheWAS trait category", y = NULL) + theme_ipf_lf() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
          axis.text.y = element_text(size = 6.5))
} else {
  pd <- ggplot() + annotate("text", x = 0, y = 0, label = "No PheWAS hits") + theme_ipf_lf()
}

## ---------- panel e: power / MDE ----------
pe <- ggplot(m[p < 0.1], aes(x = abs(mde80), y = abs(b), color = layer_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_color_manual(values = layer_cols, name = "Evidence layer") +
  labs(x = "MDE (80% power)", y = "|MR beta|") + theme_ipf_lf()

fig4 <- (pa | pb) / (pc | pd | pe) + plot_annotation(tag_levels = "a")
save_panel(fig4, file.path(OUT, "fig4_main"), width = 14, height = 10)

## ---------- supplements ----------
layer_short <- c(GTEx_liver_eQTL = "GTEx liver", eQTLGen_blood_eQTL = "eQTLGen blood",
                 INTERVAL_pQTL = "INTERVAL pQTL")
instc[, lg_short := layer_short[layer_group]]
s1a <- ggplot(instc, aes(x = lg_short, y = f_stat, fill = layer_group)) +
  geom_violin(scale = "width", color = NA, alpha = 0.75) +
  geom_boxplot(width = 0.15, outlier.size = 0.4) +
  scale_y_log10() +
  scale_fill_manual(values = layer_cols, guide = "none") +
  labs(x = NULL, y = "Instrument F statistic (log10)") + theme_ipf_lf()
st <- m[!is.na(steiger_ok), .(n = .N, correct = sum(steiger_ok, na.rm = TRUE)), by = layer_group]
st[, `:=`(frac = correct / n, lg_short = layer_short[layer_group])]
s1b <- ggplot(st, aes(x = lg_short, y = frac)) + geom_col(fill = "#BEF3C5", alpha = 0.75, width = 0.6) +
  geom_text(aes(label = paste0(correct, "/", n)), vjust = -0.3, size = 3, family = "Liberation Sans") +
  ylim(0, 1.15) + labs(x = NULL, y = "Steiger correct-direction fraction") + theme_ipf_lf()
save_panel(s1a | s1b, file.path(OUT, "fig4S1_instrument_qc"), width = 9, height = 4)

eg <- m[!is.na(egger_int_p)]
if (nrow(eg)) {
  s2 <- ggplot(eg, aes(x = egger_int, y = -log10(p), color = layer_group)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
    geom_point(size = 1.6, alpha = 0.75) + scale_color_manual(values = layer_cols, name = "Layer") +
    labs(x = "MR-Egger intercept", y = "MR -log10(p)") + theme_ipf_lf()
  save_panel(s2, file.path(OUT, "fig4S2_egger"), width = 6, height = 4.5)
}

if (nrow(ann)) {
  al <- melt(ann, id.vars = "gene", measure.vars = patterns("^sceqtl_"),
             variable.name = "ct", value.name = "hit")
  al[, ct := sub("sceqtl_", "", ct)]
  al <- al[gene %in% ann[any_sceqtl == TRUE]$gene]
  if (nrow(al)) {
    s3 <- ggplot(al, aes(x = ct, y = gene, fill = hit)) +
      geom_tile(color = "white") +
      scale_fill_manual(values = c(`TRUE` = "#BFDA8E", `FALSE` = "#FAF9F3"), name = "sc-eQTL") +
      labs(x = "Liver cell type (GSE289173)", y = NULL) + theme_ipf_lf() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 6))
    save_panel(s3, file.path(OUT, "fig4S3_sceqtl"), width = 5.5,
               height = max(3, 0.25 * length(unique(al$gene)) + 1))
  }
}

pcc <- fread(file.path(QC, "positive_control_results.csv"))
if (nrow(pcc) && "b" %in% names(pcc)) {
  pcc[, lab := paste(gene, "|", outcome)]
  pcc <- pcc[order(b)][, lab := factor(lab, levels = unique(lab))]
  s4 <- ggplot(pcc, aes(x = b, y = lab)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.3) +
    geom_point(aes(color = gene), size = 2, alpha = 0.75) +
    scale_color_manual(values = c(HSD17B13 = "#DCA7EB", MBOAT7 = "#BFDA8E"), guide = "none") +
    labs(x = "MR effect (eQTLGen blood eQTL)", y = NULL) + theme_ipf_lf() +
    theme(axis.text.y = element_text(size = 6))
  save_panel(s4, file.path(OUT, "fig4S4_positive_control"), width = 7, height = 4)
}
cat("FIGS_DONE\n")
