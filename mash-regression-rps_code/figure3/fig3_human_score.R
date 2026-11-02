#!/usr/bin/env Rscript
# Human: component UCell scoring + AUCell robustness + MASLD vs Healthy tests
suppressMessages({ library(Seurat); library(dplyr); library(UCell); library(AUCell) })
OUT <- "/workspace/fig3_out"
h <- readRDS(file.path(OUT, "human_sc_integrated.rds"))
h <- JoinLayers(h)

rev_g <- intersect(readLines(file.path(OUT, "genes_reversal_human.txt")), rownames(h))
eme_g <- intersect(readLines(file.path(OUT, "genes_emergent_human.txt")), rownames(h))
cat("human genes: reversal", length(rev_g), "emergent", length(eme_g), "\n")
d <- GetAssayData(h, layer = "data")
sc <- as.data.frame(ScoreSignatures_UCell(d, features = list(reversal = rev_g, emergent = eme_g), name = "_prog"))
stopifnot(identical(rownames(sc), colnames(h)))
h@meta.data <- cbind(h@meta.data, sc)

set.seed(1)
cells <- colnames(d); if (length(cells) > 20000) cells <- sample(cells, 20000)
rankings <- AUCell_buildRankings(d[, cells], plotStats = FALSE, verbose = FALSE)
auc <- AUCell_calcAUC(list(reversal = rev_g, emergent = eme_g), rankings,
                      aucMaxRank = ceiling(0.05 * nrow(rankings)), verbose = FALSE)
auc_mat <- getAUC(auc)
write.csv(data.frame(cell = colnames(auc_mat), reversal_AUCell = auc_mat["reversal", ],
                     emergent_AUCell = auc_mat["emergent", ]),
          file.path(OUT, "human_aucell_scores.csv"), row.names = FALSE)

saveRDS(h, file.path(OUT, "human_sc_scored_tmp.rds"))
file.rename(file.path(OUT, "human_sc_scored_tmp.rds"), file.path(OUT, "human_sc_scored.rds"))
write.csv(h@meta.data, file.path(OUT, "human_sc_metadata.csv"))

md <- h@meta.data
summ <- bind_rows(lapply(c("reversal", "emergent"), function(comp) {
  data.frame(celltype = md$celltype_label, grp = md$group, comp = comp,
             score = md[[paste0(comp, "_prog")]])
})) %>% group_by(comp, celltype, grp) %>%
  summarise(median = median(score), mean = mean(score), n = n(), .groups = "drop")
write.csv(summ, file.path(OUT, "ucell_summary_human.csv"), row.names = FALSE)

tests <- bind_rows(lapply(c("reversal", "emergent"), function(comp) {
  v <- paste0(comp, "_prog")
  bind_rows(lapply(unique(md$celltype_label), function(ct) {
    sub <- md[md$celltype_label == ct, ]
    if (sum(sub$group == "Healthy") < 10 || sum(sub$group == "MASLD") < 10)
      return(data.frame(comp = comp, celltype = ct, p = NA_real_, med_healthy = NA_real_, med_masld = NA_real_))
    data.frame(comp = comp, celltype = ct,
               p = wilcox.test(sub[[v]] ~ sub$group)$p.value,
               med_healthy = median(sub[[v]][sub$group == "Healthy"]),
               med_masld = median(sub[[v]][sub$group == "MASLD"]))
  }))
}))
tests$fdr <- p.adjust(tests$p, "BH")
write.csv(tests, file.path(OUT, "ucell_tests_human.csv"), row.names = FALSE)
tests <- tests[!is.na(tests$p), ]; tests <- tests[order(tests$fdr), ]
writeLines(sprintf("%-9s %-22s medHealthy=%.3f medMASLD=%.3f p=%.1e fdr=%.1e",
                   tests$comp, tests$celltype, tests$med_healthy, tests$med_masld, tests$p, tests$fdr))
cat("DONE\n")
