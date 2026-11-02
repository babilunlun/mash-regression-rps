#!/usr/bin/env Rscript
# make_fig5.R — Figure 5 main (6 panels) + supplements: tractability of causal targets
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table); library(ggplot2); library(patchwork)
source("/workspace/fig1_theme.R")
OUT <- "/workspace/fig5/figs"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
QC4 <- "/workspace/fig4/qc"

tt <- fread("/workspace/fig5/fig5_target_table.csv")
drugs <- fread("/workspace/fig5/known_drugs.csv")
meta <- fread(file.path(QC4, "evidence_tiers.csv"))
sus <- fread(file.path(QC4, "susie_coloc_h3_recheck.csv"))
m4 <- fread(file.path(QC4, "mr_results_annotated.csv"))

leads <- c("TP53INP1", "SLAMF6")
scorecard8 <- c("TP53INP1", "SLAMF6", "FN3K", "MYO18A", "OAS2", "TSPAN3", "CTSS", "ICAM1")

## genetic evidence score (core quartet, mirrors Fig4 panel f: best row per gene by tier then FDR)
sus_ok <- sus[verdict == "colocalised_signal_pair", unique(gene)]
setorder(meta, tier, best_fdr)
gs <- meta[, .SD[1], by = gene][, .(gene, clin_p, PPH4_max, modality, best_fdr, tier)]
gs[, `:=`(coloc_ok = (!is.na(PPH4_max) & PPH4_max > 0.8) | gene %in% sus_ok,
          clin_ok = !is.na(clin_p) & clin_p < 0.05,
          dir_ok = modality %in% c("agonist", "inhibitor"),
          fdr_ok = !is.na(best_fdr) & best_fdr < 0.05)]
gs[, gen_score := fdr_ok + coloc_ok + clin_ok + dir_ok]
tt <- merge(tt, gs[, .(gene, gen_score, tier, modality)], by = "gene", all.x = TRUE)
tt[, tract_frac := pmax(SM_n_true / SM_n_total, AB_n_true / AB_n_total,
                        PR_n_true / PR_n_total, OC_n_true / OC_n_total, na.rm = TRUE)]
tt[, cls := fcase(gene %in% leads, "Lead",
                  gene %in% scorecard8, "Scorecard",
                  default = "FDR hit")]
cls_cols <- c(Lead = "#E2EA4D", Scorecard = "#A7C0DF", `FDR hit` = "#C8CACB")

## ---------- panel a: genetics x tractability matrix ----------
pa <- ggplot(tt, aes(x = gen_score, y = tract_frac)) +
  geom_jitter(aes(color = cls, size = pmin(n_known_drugs, 4)), width = 0.18, height = 0.02, alpha = 0.75) +
  ggrepel::geom_text_repel(data = tt[cls != "FDR hit" | n_known_drugs > 0],
                           aes(label = gene), size = 2.6, family = "Liberation Sans",
                           show.legend = FALSE, max.overlaps = 25, seed = 7) +
  scale_color_manual(values = cls_cols, name = NULL) +
  scale_size_continuous(name = "Known drugs", range = c(1.5, 5)) +
  scale_x_continuous(breaks = 0:4) +
  labs(x = "Genetic evidence score (0-4 core criteria)", y = "Best tractability fraction (Open Targets)") +
  theme_ipf_lf()

## ---------- panel b: tractability heatmap, scorecard 8 ----------
th <- melt(tt[gene %in% scorecard8],
           id.vars = "gene", measure.vars = c("SM_n_true", "AB_n_true", "PR_n_true", "OC_n_true"),
           variable.name = "modality", value.name = "n")
tot <- melt(tt[gene %in% scorecard8],
            id.vars = "gene", measure.vars = c("SM_n_total", "AB_n_total", "PR_n_total", "OC_n_total"),
            variable.name = "modality2", value.name = "nt")
th[, frac := n / tot$nt]
th[, modality := factor(sub("_n_true", "", modality), levels = c("SM", "AB", "PR", "OC"),
                        labels = c("Small molecule", "Antibody", "PROTAC", "Other"))]
gord <- tt[gene %in% scorecard8][order(-gen_score, -tract_frac)]$gene
th[, gene := factor(gene, levels = rev(gord))]
pb <- ggplot(th, aes(x = modality, y = gene, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = n), size = 3, family = "Liberation Sans") +
  scale_fill_gradient(low = "#FAF9F3", high = "#BFDA8E", name = "Bucket\nfraction") +
  labs(x = "Modality", y = NULL) + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5))

## ---------- panel c: safety — LOEUF vs tissue breadth ----------
lab_c <- tt[cls != "FDR hit" | oe_lof_upper < 0.35]
pc <- ggplot(tt[!is.na(oe_lof_upper)], aes(x = oe_lof_upper, y = n_tissues_tpm1)) +
  geom_vline(xintercept = 0.35, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_point(aes(color = cls), size = 2, alpha = 0.75) +
  ggrepel::geom_text_repel(data = lab_c[!is.na(oe_lof_upper)],
                           aes(label = gene), color = "grey25", size = 2.6, family = "Liberation Sans",
                           show.legend = FALSE, max.overlaps = 30, seed = 7,
                           force_pull = 0.5, min.segment.length = 0.2) +
  scale_color_manual(values = cls_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.03))) +
  coord_cartesian(ylim = c(NA, 58)) +
  labs(x = "LOEUF (LoF intolerance; left = constrained)", y = "Tissues with median TPM >= 1") +
  theme_ipf_lf()

## ---------- panel d: human liver cell-type expression (MASLD) ----------
sc <- fread("/workspace/fig5/sc_expr_human_8genes.csv")
scm <- sc[group == "MASLD", .(mean_expr = mean(mean_expr)), by = .(gene, lineage)]
scm[, gene := factor(gene, levels = rev(gord))]
pd <- ggplot(scm, aes(x = lineage, y = gene, fill = mean_expr)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = "#FAF9F3", high = "#DCA7EB", name = "Mean log-norm\nexpression", trans = "sqrt") +
  labs(x = "Liver cell lineage (MASLD)", y = NULL) + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

## ---------- panel e: known drugs ----------
dg <- merge(drugs, tt[, .(gene, cls, gen_score, modality)], by = "gene")
# verified corrections (2026-10) to Open Targets maxClinicalStage:
# sabatolimab — Ph3 STIMULUS-MDS2 failed OS endpoint, program discontinued (never approved)
# glembatumumab vedotin — Ph2b METRIC failed 2018, program discontinued (never approved)
dg[drug == "SABATOLIMAB", stage := "PHASE_3"]
dg[drug == "GLEMBATUMUMAB VEDOTIN", stage := "PHASE_2"]
dg[, stage := factor(stage, levels = c("APPROVAL", "PHASE_3", "PHASE_2_3", "PHASE_2", "IND", "UNKNOWN"))]
dg[, drug_tc := tools::toTitleCase(tolower(drug))]
# one row per drug: labels on the axis, no in-tile text -> no clipping/overlap
dg[, ylab := paste0(drug_tc, " (", gene, ")")]
dg[gene == "ATP1B3", ylab := paste0(ylab, " \u2021")]  # dagger: complex-level annotation, canonical target is ATP1A1
dg[, dir_ok := modality %in% c("agonist", "inhibitor")]
dg[, dir_lab := fifelse(dir_ok, substr(modality, 1, 3), "n.d.")]
dg[, cls_rank := fifelse(cls == "Scorecard", 0, 1)]
setorder(dg, cls_rank, -gen_score, gene, stage, drug_tc)
dg[, ylab := factor(ylab, levels = rev(unique(ylab)))]
pe <- ggplot(dg, aes(x = stage, y = ylab, fill = cls)) +
  geom_tile(color = "white", linewidth = 0.4, alpha = 0.75, height = 0.8) +
  scale_fill_manual(values = cls_cols, name = NULL) +
  scale_x_discrete(drop = FALSE, labels = c("APPROVAL" = "Approved", "PHASE_3" = "Ph3", "PHASE_2_3" = "Ph2/3",
                                            "PHASE_2" = "Ph2", "IND" = "IND", "UNKNOWN" = "Unk.")) +
  labs(x = "Max clinical stage (any indication)", y = NULL) + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 5.5))
# direction-of-effect strip: does the drug's MoA match the genetically predicted direction?
pdir <- ggplot(dg, aes(x = "x", y = ylab, fill = dir_ok)) +
  geom_tile(color = "white", linewidth = 0.4, alpha = 0.75, height = 0.8) +
  geom_text(aes(label = dir_lab), size = 1.6, family = "Liberation Sans") +
  scale_fill_manual(values = c(`TRUE` = "#BFDA8E", `FALSE` = "#C8CACB"), guide = "none") +
  scale_x_discrete(labels = c("x" = "Genetic\ndirection")) +
  labs(x = NULL, y = NULL) + theme_ipf_lf() +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6))
# wrap as one grob so the composite keeps a single panel tag
pe <- wrap_elements(full = pe + pdir + plot_layout(widths = c(8, 1)))

## ---------- panel f: modality scorecard for the 2 leads ----------
sc_json <- function(g) {
  s <- jsonlite::fromJSON(paste0("/workspace/fig5/", g, "_scorecard.json"))
  data.table(gene = g, modality = names(s$modalities),
             overall = sapply(s$modalities, function(x) x$overall))
}
mf <- rbindlist(lapply(leads, sc_json))
mf[, modality := factor(modality, levels = c("SM", "AB", "PR", "OC"),
                        labels = c("Small molecule", "Antibody", "PROTAC", "Other"))]
pf <- ggplot(mf, aes(x = modality, y = gene, fill = overall)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = overall), size = 3.2, family = "Liberation Sans") +
  scale_fill_gradient(low = "#FAF9F3", high = "#A7C0DF", limits = c(0, 3), name = "Modality\nscore (0-3)") +
  labs(x = "Modality (protein-targeting)", y = NULL) + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5))

fig5 <- (pa | pb) / (pc | pd) / (pe | pf) + plot_annotation(tag_levels = "a")
save_panel(fig5, file.path(OUT, "fig5_main"), width = 13, height = 13.5)

## ---------- S1: full 48-gene tractability heatmap ----------
th48 <- melt(tt, id.vars = c("gene", "gen_score"),
             measure.vars = c("SM_n_true", "AB_n_true", "PR_n_true", "OC_n_true"),
             variable.name = "modality", value.name = "n")
tot48 <- melt(tt, id.vars = "gene",
              measure.vars = c("SM_n_total", "AB_n_total", "PR_n_total", "OC_n_total"),
              value.name = "nt")
th48[, frac := n / tot48$nt]
th48[, modality := factor(sub("_n_true", "", modality), levels = c("SM", "AB", "PR", "OC"),
                          labels = c("Small molecule", "Antibody", "PROTAC", "Other"))]
gord48 <- tt[order(-gen_score, -tract_frac)]$gene
th48[, gene := factor(gene, levels = rev(gord48))]
s1 <- ggplot(th48, aes(x = modality, y = gene, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#FAF9F3", high = "#BFDA8E", name = "Bucket\nfraction") +
  labs(x = "Modality", y = NULL) + theme_ipf_lf() +
  theme(axis.text.y = element_text(size = 5.5), axis.text.x = element_text(angle = 30, hjust = 1))
save_panel(s1, file.path(OUT, "fig5S1_tractability_48"), width = 6, height = 9)

## ---------- S2: GTEx tissue heatmap for scorecard 8 ----------
gt <- fread("/workspace/fig5/gtex_tissue_tpm.csv")
gt8 <- gt[gene %in% scorecard8]
gt8[, gene := factor(gene, levels = rev(gord))]
gt8[, lt := log10(median_tpm + 0.1)]
s2 <- ggplot(gt8, aes(x = tissue, y = gene, fill = lt)) +
  geom_tile(color = NA) +
  scale_fill_gradient(low = "#FAF9F3", high = "#A7C0C7", name = "log10(TPM+0.1)") +
  labs(x = "GTEx tissue", y = NULL) + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 5))
save_panel(s2, file.path(OUT, "fig5S2_gtex_tissues"), width = 10, height = 4)

## ---------- S3: mouse sc expression ----------
scem <- fread("/workspace/fig5/sc_expr_mouse_8genes.csv")
scem[, gene := tools::toTitleCase(tolower(gene))]
scem$gene[scem$gene == "Trp53inp1"] <- "TP53INP1"; scem$gene[scem$gene == "Slamf6"] <- "SLAMF6"
scem$gene[scem$gene == "Fn3k"] <- "FN3K"; scem$gene[scem$gene == "Myo18a"] <- "MYO18A"
scem$gene[scem$gene == "Oas2"] <- "OAS2"; scem$gene[scem$gene == "Tspan3"] <- "TSPAN3"
scem$gene[scem$gene == "Ctss"] <- "CTSS"; scem$gene[scem$gene == "Icam1"] <- "ICAM1"
scemm <- scem[, .(mean_expr = mean(mean_expr)), by = .(gene, lineage, group)]
# complete grid: lineages absent in a condition (e.g. HSC_MF in FATMASH_reg2wk) render as NA tiles, not gaps
scemm <- merge(CJ(gene = unique(scemm$gene), lineage = unique(scemm$lineage), group = unique(scemm$group)),
               scemm, by = c("gene", "lineage", "group"), all.x = TRUE)
scemm[, gene := factor(gene, levels = rev(gord))]
s3 <- ggplot(scemm, aes(x = lineage, y = gene, fill = mean_expr)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~group, ncol = 3) +
  scale_fill_gradient(low = "#FAF9F3", high = "#DCA7EB", name = "Mean log-norm\nexpression", na.value = "grey92") +
  labs(x = "Mouse liver cell lineage", y = NULL) + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 6.5))
save_panel(s3, file.path(OUT, "fig5S3_mouse_sc"), width = 11, height = 6)

## ---------- S4: LOEUF distribution ----------
s4 <- ggplot(tt[!is.na(oe_lof_upper)], aes(x = reorder(gene, oe_lof_upper), y = oe_lof_upper, fill = cls)) +
  geom_hline(yintercept = 0.35, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_col(alpha = 0.75, width = 0.7) +
  scale_fill_manual(values = cls_cols, name = NULL) +
  coord_flip() + labs(x = NULL, y = "LOEUF (gnomAD v4)") + theme_ipf_lf() +
  theme(axis.text.y = element_text(size = 5))
save_panel(s4, file.path(OUT, "fig5S4_loeuf"), width = 6, height = 9)

fwrite(tt, "/workspace/fig5/fig5_target_table_full.csv")
dg_out <- dg[, .(gene, drug, type, stage = as.character(stage), moa, genetic_direction = modality)]
dg_out[, direction_match := fcase(genetic_direction %in% c("agonist", "inhibitor"), "yes",
                                  default = "direction unresolved")]
dg_out[drug == "SABATOLIMAB", note := "OT maxClinicalStage=APPROVAL incorrect; Ph3 STIMULUS-MDS2 failed OS endpoint, development discontinued (never approved)"]
dg_out[drug == "GLEMBATUMUMAB VEDOTIN", note := "OT maxClinicalStage=APPROVAL incorrect; Ph2b METRIC failed 2018, development discontinued (never approved)"]
dg_out[gene == "ATP1B3", note := "Cardiac glycosides' canonical target is the catalytic alpha subunit ATP1A1; ATP1B3 (beta-3 regulatory subunit) link is likely complex-level annotation. Narrow therapeutic index - not a realistic repurposing candidate."]
fwrite(dg_out, "/workspace/fig5/known_drugs_corrected.csv")
cat("FIG5_DONE\n")
