#!/usr/bin/env Rscript
# Fig3-S9: human SAMac <-> reversal component circularity sensitivity (leave-markers-out)
# Mirrors Fig3-S8 (mouse). Predictor = SAMac_OLR1 7-gene score (unchanged);
# outcome = reversal component WITHOUT any SAMac marker gene (65 -> 61).
suppressMessages({ library(Seurat); library(dplyr); library(UCell) })
OUT <- "/workspace/fig3_out"

h <- readRDS(file.path(OUT, "human_sc_scored.rds"))
h <- JoinLayers(h)
md <- h@meta.data
cat("total cells:", nrow(md), "(expect 70424)\n")

## ---- gene sets ----
rev_full <- intersect(readLines(file.path(OUT, "genes_reversal_human.txt")), rownames(h))
samac_score_genes <- intersect(c("OLR1", "TREM2", "CD9", "SPP1", "GPNMB", "LGALS3", "FABP5"),
                               rownames(h))
samac_all <- intersect(unique(c(samac_score_genes, "CSTB", "ACP5", "CTSF")), rownames(h))
overlap <- intersect(rev_full, samac_all)
rev_lo <- setdiff(rev_full, samac_all)
cat("reversal full:", length(rev_full), "| SAMac markers in component:", length(overlap),
    "(", paste(overlap, collapse = ", "), ") | leave-out n:", length(rev_lo), "\n")

## ---- re-score ----
d <- GetAssayData(h, layer = "data")
sc <- as.data.frame(ScoreSignatures_UCell(
  d, features = list(reversalLO = rev_lo, SAMac = samac_score_genes), name = "_x"))
stopifnot(identical(rownames(sc), rownames(md)))
md$reversalLO_prog <- sc$reversalLO_x
md$SAMac_sc <- sc$SAMac_x

## ---- test A: cell-type ranking full vs leave-out, per group (Fig3b human claim) ----
rankA <- md %>% group_by(group, celltype_label) %>%
  summarise(median_full = median(reversal_prog), median_lo = median(reversalLO_prog),
            n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(rank_full = rank(-median_full), rank_lo = rank(-median_lo)) %>%
  arrange(group, rank_lo)
print(as.data.frame(rankA[rankA$group == "MASLD", ]))
write.csv(rankA, file.path(OUT, "fig3s9_celltype_ranking_leaveout_human.csv"), row.names = FALSE)

## ---- test B: per-group Spearman cor(SAMac score, reversal) within myeloid, full vs LO ----
mye_types <- c("Kupffer", "SAMac (OLR1+)", "Monocyte/MF", "Neutrophil", "cDC")
mye <- md[md$celltype_label %in% mye_types, ]
cat("human myeloid cells:", nrow(mye), "\n")
corB <- mye %>% group_by(group) %>%
  summarise(n = n(),
            rho_full = cor(SAMac_sc, reversal_prog, method = "spearman"),
            rho_lo   = cor(SAMac_sc, reversalLO_prog, method = "spearman"),
            p_full   = suppressWarnings(cor.test(SAMac_sc, reversal_prog, method = "spearman"))$p.value,
            p_lo     = suppressWarnings(cor.test(SAMac_sc, reversalLO_prog, method = "spearman"))$p.value,
            .groups = "drop")
print(as.data.frame(corB))
write.csv(corB, file.path(OUT, "fig3s9_cor_full_vs_leaveout_human.csv"), row.names = FALSE)

## ---- test C: SAMac-high quartile vs rest (leave-out), per group ----
hiC <- mye %>% group_by(group) %>%
  mutate(samac_hi = SAMac_sc >= quantile(SAMac_sc, 0.75)) %>%
  summarise(p_lo = wilcox.test(reversalLO_prog ~ samac_hi)$p.value,
            med_hi = median(reversalLO_prog[samac_hi]),
            med_lo = median(reversalLO_prog[!samac_hi]), .groups = "drop")
print(as.data.frame(hiC))
write.csv(hiC, file.path(OUT, "fig3s9_samac_hi_tests_leaveout_human.csv"), row.names = FALSE)

## ---- test D: SAMac cluster MASLD vs Healthy, full vs leave-out ----
sam <- md[md$celltype_label == "SAMac (OLR1+)", ]
testD <- data.frame(
  n_healthy = sum(sam$group == "Healthy"), n_masld = sum(sam$group == "MASLD"),
  full_healthy = median(sam$reversal_prog[sam$group == "Healthy"]),
  full_masld = median(sam$reversal_prog[sam$group == "MASLD"]),
  p_full = wilcox.test(sam$reversal_prog ~ sam$group)$p.value,
  lo_healthy = median(sam$reversalLO_prog[sam$group == "Healthy"]),
  lo_masld = median(sam$reversalLO_prog[sam$group == "MASLD"]),
  p_lo = wilcox.test(sam$reversalLO_prog ~ sam$group)$p.value)
print(testD)
write.csv(testD, file.path(OUT, "fig3s9_samac_masld_vs_healthy_human.csv"), row.names = FALSE)

saveRDS(md, file.path(OUT, "human_metadata_leaveout.rds"))
cat("DONE\n")
