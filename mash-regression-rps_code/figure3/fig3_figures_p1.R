#!/usr/bin/env Rscript
# Figure 3 phase-1 figures: panel a (human UMAP), panel b (UCell scores human+mouse),
# panel c (mouse dynamics), S1 (QC), S2 (UCell vs AUCell)
suppressMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(tidyr) })
source("/workspace/fig1_theme.R")   # theme_ipf_lf, palette8, save_panel
OUT <- "/workspace/fig3_out"; FIG <- "/workspace/fig3_fig"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

h <- readRDS(file.path(OUT, "human_sc_scored.rds"))
mo <- readRDS(file.path(OUT, "mouse_sc_scored.rds"))

ct_levels <- c("Hepatocyte", "Cholangiocyte", "Endothelial/LSEC", "HSC/myofibroblast",
               "Kupffer", "SAMac (OLR1+)", "SAMac (Trem2+)", "Monocyte/MF", "Neutrophil",
               "cDC", "pDC", "Mast cell", "T/NK", "B cell", "Plasma", "Platelet",
               "Erythrocyte", "Other")
# palette8 for the 8 core lineages; muted extensions for the rest (documented in README)
ct_pal <- c(Hepatocyte = "#DCA7EB", Cholangiocyte = "#B1C4D7", `Endothelial/LSEC` = "#BFDA8E",
            `HSC/myofibroblast` = "#E2EA4D", Kupffer = "#BEF3C5", `SAMac (OLR1+)` = "#C7B0C7",
            `SAMac (Trem2+)` = "#C7B0C7", `Monocyte/MF` = "#A7C0DF", Neutrophil = "#FDB863",
            cDC = "#F4A582", pDC = "#80CDC1", `Mast cell` = "#B2ABD2", `T/NK` = "#C8CACB",
            `B cell` = "#E78AC3", Plasma = "#999999", Platelet = "#DDDDDD",
            Erythrocyte = "#F0F0F0", Other = "#F0F0F0")

## ---------- panel a: human UMAP ----------
pa <- DimPlot(h, group.by = "celltype_label", cols = ct_pal, label = TRUE, label.size = 2.6,
              repel = TRUE, pt.size = 0.05) +
  theme_ipf_lf() + ggtitle("Human liver scRNA (GSE298719)") +
  theme(legend.position = "none")
save_panel(pa, file.path(FIG, "fig3a_human_umap"), width = 6.5, height = 5)

## ---------- panel b: UCell per cell type ----------
mk_long <- function(obj, group_var, species) {
  md <- obj@meta.data
  bind_rows(lapply(c("reversal", "emergent"), function(comp) {
    data.frame(species = species, comp = comp, celltype = md$celltype_label,
               grp = md[[group_var]], score = md[[paste0(comp, "_prog")]])
  }))
}
dl <- bind_rows(mk_long(h, "group", "Human"), mk_long(mo, "grp2", "Mouse"))
dl <- dl %>% filter(!celltype %in% c("Other", "Erythrocyte", "Platelet"))
dl$grp <- factor(dl$grp, levels = c("Healthy", "MASLD", "peak", "reg2wk", "reg6wk"))
grp_pal <- c(Healthy = "#C8CACB", MASLD = "#E2EA4D", peak = "#E2EA4D", reg2wk = "#B1C4D7", reg6wk = "#BFDA8E")

pb <- ggplot(dl %>% filter(grp %in% c("Healthy", "MASLD", "peak", "reg6wk")),
             aes(x = celltype, y = score, fill = grp)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.2, linewidth = 0.2, position = position_dodge(0.75)) +
  facet_grid(comp ~ species, scales = "free") +
  scale_fill_manual(values = grp_pal) +
  theme_ipf_lf() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  labs(x = NULL, y = "UCell score", fill = NULL)
save_panel(pb, file.path(FIG, "fig3b_ucell_scores"), width = 11, height = 6)

## ---------- panel c: mouse dynamics ----------
prop <- read.csv(file.path(OUT, "mouse_celltype_proportions.csv"))
prop <- prop %>% filter(!celltype_label %in% c("Other", "Erythrocyte", "Platelet"))
prop$timepoint <- factor(prop$timepoint, levels = c("peak", "reg2wk", "reg6wk"))
pc1 <- ggplot(prop, aes(x = timepoint, y = frac, fill = celltype_label)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  facet_wrap(~model) + scale_fill_manual(values = ct_pal) +
  theme_ipf_lf() + labs(x = NULL, y = "Cell fraction", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
traj <- read.csv(file.path(OUT, "mouse_score_trajectories.csv"))
traj <- traj %>% filter(!celltype_label %in% c("Other", "Erythrocyte", "Platelet"))
traj$timepoint <- factor(traj$timepoint, levels = c("peak", "reg2wk", "reg6wk"))
traj_l <- traj %>% pivot_longer(c(reversal, emergent), names_to = "comp", values_to = "score")
pc2 <- ggplot(traj_l, aes(x = timepoint, y = score, group = celltype_label, color = celltype_label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1) +
  facet_grid(comp ~ model, scales = "free_y") + scale_color_manual(values = ct_pal) +
  theme_ipf_lf() + labs(x = NULL, y = "Median UCell", color = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), legend.position = "none")
pc <- pc1 + pc2 + plot_layout(widths = c(1, 1.4))
save_panel(pc, file.path(FIG, "fig3c_mouse_dynamics"), width = 11, height = 5)

## ---------- S1: QC ----------
qc_h <- read.csv(file.path(OUT, "qc_human_sc_raw.csv")); qc_h$species <- "Human"
qc_m <- read.csv(file.path(OUT, "qc_mouse_sc_raw.csv")); qc_m$species <- "Mouse"
qc <- bind_rows(qc_h, qc_m)
qc_l <- qc %>% pivot_longer(c(median_nFeature, median_nCount, median_mt), names_to = "metric", values_to = "value")
s1a <- ggplot(qc_l, aes(x = sample, y = value, fill = species)) +
  geom_col(width = 0.7) + facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(Human = palette8[1], Mouse = palette8[2])) +
  theme_ipf_lf() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4), legend.position = "top") +
  labs(x = NULL, y = NULL)
s1b <- DimPlot(h, group.by = "sample", pt.size = 0.05) + theme_ipf_lf() +
  ggtitle("Human: by sample") + theme(legend.position = "none")
s1c <- DimPlot(mo, group.by = "condition", pt.size = 0.05) + theme_ipf_lf() +
  ggtitle("Mouse: by condition") + theme(legend.position = "none")
s1 <- s1a | (s1b / s1c)
save_panel(s1, file.path(FIG, "fig3S1_qc"), width = 11, height = 7)

## ---------- S2: UCell vs AUCell ----------
mk_auc <- function(auc_csv, obj, species) {
  a <- read.csv(auc_csv); md <- obj@meta.data
  a$celltype <- md$celltype_label[match(a$cell, rownames(md))]
  a <- a %>% filter(!is.na(celltype))
  u <- data.frame(cell = rownames(md), reversal_prog = md$reversal_prog,
                  emergent_prog = md$emergent_prog, celltype = md$celltype_label)
  j <- inner_join(a, u, by = c("cell", "celltype"))
  list(species = species,
       cor_rev = cor(j$reversal_AUCell, j$reversal_prog, method = "spearman"),
       cor_eme = cor(j$emergent_AUCell, j$emergent_prog, method = "spearman"),
       plot = j)
}
r_h <- mk_auc(file.path(OUT, "human_aucell_scores.csv"), h, "Human")
r_m <- mk_auc(file.path(OUT, "mouse_aucell_scores.csv"), mo, "Mouse")
cat("UCell-AUCell Spearman: human rev", r_h$cor_rev, "eme", r_h$cor_eme,
    "| mouse rev", r_m$cor_rev, "eme", r_m$cor_eme, "\n")
j_all <- bind_rows(r_h$plot %>% mutate(species = "Human"), r_m$plot %>% mutate(species = "Mouse"))
s2a <- ggplot(j_all, aes(x = reversal_prog, y = reversal_AUCell)) +
  geom_point(size = 0.05, alpha = 0.15, color = palette8[1]) +
  facet_wrap(~species, scales = "free") +
  theme_ipf_lf() + scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  labs(x = "UCell (reversal)", y = "AUCell (reversal)")
s2b <- ggplot(j_all, aes(x = emergent_prog, y = emergent_AUCell)) +
  geom_point(size = 0.05, alpha = 0.15, color = palette8[3]) +
  facet_wrap(~species, scales = "free") +
  theme_ipf_lf() + scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  labs(x = "UCell (emergent)", y = "AUCell (emergent)")
s2 <- s2a | s2b
save_panel(s2, file.path(FIG, "fig3S2_ucell_vs_aucell"), width = 9, height = 3.8)
cat("FIGURES DONE\n")
