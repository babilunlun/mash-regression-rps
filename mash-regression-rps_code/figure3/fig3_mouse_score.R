#!/usr/bin/env Rscript
# Mouse-only: component UCell scoring + AUCell robustness + dynamics tables
suppressMessages({ library(Seurat); library(dplyr); library(UCell); library(AUCell) })
OUT <- "/workspace/fig3_out"
mo <- readRDS(file.path(OUT, "mouse_sc_integrated.rds"))
mo <- JoinLayers(mo)

rev_g <- intersect(readLines(file.path(OUT, "genes_reversal_mouse.txt")), rownames(mo))
eme_g <- intersect(readLines(file.path(OUT, "genes_emergent_mouse.txt")), rownames(mo))
cat("mouse genes: reversal", length(rev_g), "emergent", length(eme_g), "\n")
d <- GetAssayData(mo, layer = "data")
sc <- as.data.frame(ScoreSignatures_UCell(d, features = list(reversal = rev_g, emergent = eme_g), name = "_prog"))
stopifnot(identical(rownames(sc), colnames(mo)))
mo@meta.data <- cbind(mo@meta.data, sc)

# AUCell robustness (20k cell subsample)
set.seed(1)
cells <- colnames(d); if (length(cells) > 20000) cells <- sample(cells, 20000)
rankings <- AUCell_buildRankings(d[, cells], plotStats = FALSE, verbose = FALSE)
auc <- AUCell_calcAUC(list(reversal = rev_g, emergent = eme_g), rankings,
                      aucMaxRank = ceiling(0.05 * nrow(rankings)), verbose = FALSE)
auc_mat <- getAUC(auc)
auc_df <- data.frame(cell = colnames(auc_mat),
                     reversal_AUCell = auc_mat["reversal", ],
                     emergent_AUCell = auc_mat["emergent", ])
write.csv(auc_df, file.path(OUT, "mouse_aucell_scores.csv"), row.names = FALSE)

mo$grp2 <- factor(mo$timepoint, levels = c("peak", "reg2wk", "reg6wk"))
saveRDS(mo, file.path(OUT, "mouse_sc_scored.rds"))
write.csv(mo@meta.data, file.path(OUT, "mouse_sc_metadata.csv"))

# per-celltype summaries + tests (peak vs reg6wk, per model)
md <- mo@meta.data
long <- bind_rows(lapply(c("reversal", "emergent"), function(comp) {
  data.frame(celltype = md$celltype_label, model = md$model, timepoint = md$timepoint,
             comp = comp, score = md[[paste0(comp, "_prog")]])
}))
summ <- long %>% group_by(comp, celltype, model, timepoint) %>%
  summarise(median = median(score), mean = mean(score), n = n(), .groups = "drop")
write.csv(summ, file.path(OUT, "ucell_summary_mouse.csv"), row.names = FALSE)

tests <- bind_rows(lapply(unique(md$model), function(mdl) {
  sub <- md[md$model == mdl & md$timepoint %in% c("peak", "reg6wk"), ]
  bind_rows(lapply(c("reversal", "emergent"), function(comp) {
    v <- paste0(comp, "_prog")
    bind_rows(lapply(unique(sub$celltype_label), function(ct) {
      s2 <- sub[sub$celltype_label == ct, ]
      if (sum(s2$timepoint == "peak") < 5 || sum(s2$timepoint == "reg6wk") < 5)
        return(data.frame(model = mdl, comp = comp, celltype = ct, p = NA_real_,
                          med_peak = NA_real_, med_reg6wk = NA_real_))
      data.frame(model = mdl, comp = comp, celltype = ct,
                 p = wilcox.test(s2[[v]] ~ s2$timepoint)$p.value,
                 med_peak = median(s2[[v]][s2$timepoint == "peak"]),
                 med_reg6wk = median(s2[[v]][s2$timepoint == "reg6wk"]))
    }))
  }))
}))
tests$fdr <- p.adjust(tests$p, "BH")
write.csv(tests, file.path(OUT, "ucell_tests_mouse.csv"), row.names = FALSE)
print(tests %>% filter(!is.na(p)) %>% arrange(fdr) %>% head(15))

# proportions + trajectories
prop <- md %>% count(model, timepoint, celltype_label, name = "n") %>%
  group_by(model, timepoint) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.csv(prop, file.path(OUT, "mouse_celltype_proportions.csv"), row.names = FALSE)
# immune-only proportions (digestion bias robust)
imm <- md %>% filter(!celltype_label %in% c("Hepatocyte", "Endothelial/LSEC", "Cholangiocyte", "Erythrocyte"))
prop_imm <- imm %>% count(model, timepoint, celltype_label, name = "n") %>%
  group_by(model, timepoint) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.csv(prop_imm, file.path(OUT, "mouse_immune_proportions.csv"), row.names = FALSE)

ptests <- bind_rows(lapply(unique(md$model), function(mdl) {
  sub <- md[md$model == mdl & md$timepoint %in% c("peak", "reg6wk"), ]
  bind_rows(lapply(unique(sub$celltype_label), function(ct) {
    a <- sum(sub$celltype_label == ct & sub$timepoint == "peak"); b <- sum(sub$timepoint == "peak") - a
    cc <- sum(sub$celltype_label == ct & sub$timepoint == "reg6wk"); dd <- sum(sub$timepoint == "reg6wk") - cc
    ft <- tryCatch(fisher.test(matrix(c(a, b, cc, dd), nrow = 2, byrow = TRUE)), error = function(e) NULL)
    data.frame(model = mdl, celltype = ct, n_peak = a, n_reg6wk = cc,
               frac_peak = a / (a + b), frac_reg6wk = cc / (cc + dd),
               or = if (is.null(ft)) NA_real_ else ft$estimate,
               p = if (is.null(ft)) NA_real_ else ft$p.value)
  }))
}))
ptests$fdr <- p.adjust(ptests$p, "BH")
write.csv(ptests, file.path(OUT, "mouse_proportion_tests.csv"), row.names = FALSE)
print(ptests %>% arrange(fdr) %>% head(12))

traj <- md %>% group_by(model, timepoint, celltype_label) %>%
  summarise(reversal = median(reversal_prog), emergent = median(emergent_prog), n = n(), .groups = "drop")
write.csv(traj, file.path(OUT, "mouse_score_trajectories.csv"), row.names = FALSE)
cat("DONE\n")
