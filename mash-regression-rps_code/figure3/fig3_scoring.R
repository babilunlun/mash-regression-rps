#!/usr/bin/env Rscript
# Figure 3 panel b: UCell scoring of reversal/emergent components per cell type (human + mouse)
# + AUCell robustness (S2). Panel c: mouse cell-state dynamics across regression.
suppressMessages({ library(Seurat); library(dplyr); library(UCell); library(AUCell) })
OUT <- "/workspace/fig3_out"

score_obj <- function(obj, rev_f, eme_f, label) {
  rev_g <- readLines(rev_f); eme_g <- readLines(eme_f)
  feats <- rownames(obj)
  rev_g <- intersect(rev_g, feats); eme_g <- intersect(eme_g, feats)
  cat(label, "genes found: reversal", length(rev_g), "emergent", length(eme_g), "\n")
  obj <- AddModuleScore_UCell(obj, features = list(reversal = rev_g, emergent = eme_g), name = "_prog")
  # AUCell robustness on a downsampled set (max 20000 cells for speed)
  expr <- GetAssayData(obj, layer = "data")
  set.seed(1)
  cells <- colnames(expr); if (length(cells) > 20000) cells <- sample(cells, 20000)
  rankings <- AUCell_buildRankings(expr[, cells], plotStats = FALSE, verbose = FALSE)
  auc <- AUCell_calcAUC(list(reversal = rev_g, emergent = eme_g), rankings,
                        aucMaxRank = ceiling(0.05 * nrow(rankings)), verbose = FALSE)
  auc_mat <- getAUC(auc)
  auc_df <- data.frame(cell = colnames(auc_mat),
                       reversal_AUCell = auc_mat["reversal", ],
                       emergent_AUCell = auc_mat["emergent", ])
  list(obj = obj, auc = auc_df, n_rev = length(rev_g), n_eme = length(eme_g))
}

## ---- human ----
h <- readRDS(file.path(OUT, "human_sc_integrated.rds"))
res_h <- score_obj(h, file.path(OUT, "genes_reversal_human.txt"), file.path(OUT, "genes_emergent_human.txt"), "human")
h <- res_h$obj
saveRDS(h, file.path(OUT, "human_sc_scored.rds"))
write.csv(res_h$auc, file.path(OUT, "human_aucell_scores.csv"), row.names = FALSE)

## ---- mouse ----
mo <- readRDS(file.path(OUT, "mouse_sc_integrated.rds"))
res_m <- score_obj(mo, file.path(OUT, "genes_reversal_mouse.txt"), file.path(OUT, "genes_emergent_mouse.txt"), "mouse")
mo <- res_m$obj
saveRDS(mo, file.path(OUT, "mouse_sc_scored.rds"))
write.csv(res_m$auc, file.path(OUT, "mouse_aucell_scores.csv"), row.names = FALSE)

## ---- per-celltype stats (UCell): median score per cell type x group, Wilcoxon ----
ct_stats <- function(obj, group_var, label) {
  md <- obj@meta.data
  md$grp <- md[[group_var]]
  long <- bind_rows(lapply(c("reversal", "emergent"), function(comp) {
    v <- paste0(comp, "_prog")
    data.frame(celltype = md$celltype, grp = md$grp, comp = comp, score = md[[v]])
  }))
  summ <- long %>% group_by(comp, celltype, grp) %>%
    summarise(median = median(score), mean = mean(score), n = n(), .groups = "drop")
  # Wilcoxon between groups per cell type x comp
  tests <- long %>% group_by(comp, celltype) %>%
    summarise(p = tryCatch(wilcox.test(score ~ grp)$p.value, error = function(e) NA_real_),
              med1 = median(score[grp == levels(factor(grp))[1]]),
              med2 = median(score[grp == levels(factor(grp))[2]]), .groups = "drop") %>%
    mutate(fdr = p.adjust(p, "BH"))
  write.csv(summ, file.path(OUT, paste0("ucell_summary_", label, ".csv")), row.names = FALSE)
  write.csv(tests, file.path(OUT, paste0("ucell_tests_", label, ".csv")), row.names = FALSE)
  print(tests)
  invisible(list(summ = summ, tests = tests))
}
cat("== human: MASLD vs Healthy ==\n"); st_h <- ct_stats(h, "group", "human")
cat("== mouse: peak vs reg6wk (per model) ==\n")
mo$grp2 <- ifelse(mo$timepoint == "peak", "peak", ifelse(mo$timepoint == "reg6wk", "reg6wk", "reg2wk"))
st_m <- ct_stats(mo, "grp2", "mouse")

## ---- panel c: mouse proportions + score trajectories ----
md <- mo@meta.data
prop <- md %>% count(model, timepoint, celltype, name = "n") %>%
  group_by(model, timepoint) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.csv(prop, file.path(OUT, "mouse_celltype_proportions.csv"), row.names = FALSE)
# proportion change tests: peak vs reg6wk per model (Fisher exact per cell type)
ptests <- bind_rows(lapply(unique(md$model), function(mdl) {
  sub <- md[md$model == mdl & md$timepoint %in% c("peak", "reg6wk"), ]
  bind_rows(lapply(unique(sub$celltype), function(ct) {
    a <- sum(sub$celltype == ct & sub$timepoint == "peak"); b <- sum(sub$timepoint == "peak") - a
    c <- sum(sub$celltype == ct & sub$timepoint == "reg6wk"); d <- sum(sub$timepoint == "reg6wk") - c
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
    data.frame(model = mdl, celltype = ct, frac_peak = a / (a + b), frac_reg6wk = c / (c + d),
               or = ft$estimate, p = ft$p.value)
  }))
}))
ptests$fdr <- p.adjust(ptests$p, "BH")
write.csv(ptests, file.path(OUT, "mouse_proportion_tests.csv"), row.names = FALSE)
print(ptests %>% arrange(fdr) %>% head(12))
# score trajectories per cell type x timepoint
traj <- md %>% group_by(model, timepoint, celltype) %>%
  summarise(reversal = median(reversal_prog), emergent = median(emergent_prog), n = n(), .groups = "drop")
write.csv(traj, file.path(OUT, "mouse_score_trajectories.csv"), row.names = FALSE)
cat("DONE\n")
