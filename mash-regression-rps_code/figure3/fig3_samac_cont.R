#!/usr/bin/env Rscript
# SAMac continuous-score analysis (mouse): strengthens the 75-cell cluster claim
# using all myeloid cells. cor(SAMac marker score, component score) per condition,
# + reversal/emergent in SAMac-high (top quartile) vs rest, per condition.
suppressMessages({ library(Seurat); library(dplyr); library(UCell) })
OUT <- "/workspace/fig3_out"
mo <- readRDS(file.path(OUT, "mouse_sc_scored.rds"))
mo <- JoinLayers(mo)

samac_genes <- c("Trem2", "Cd9", "Spp1", "Gpnmb", "Fabp5", "Lgals3", "Lpl")
samac_genes <- intersect(samac_genes, rownames(mo))
d <- GetAssayData(mo, layer = "data")
sc <- as.data.frame(ScoreSignatures_UCell(d, features = list(SAMac = samac_genes), name = "_sc"))
mo@meta.data <- cbind(mo@meta.data, sc)

mye_types <- c("Kupffer", "SAMac (Trem2+)", "Monocyte/MF", "Neutrophil", "cDC", "pDC")
md <- mo@meta.data[mo@meta.data$celltype_label %in% mye_types, ]
cat("myeloid cells:", nrow(md), "\n")

# (1) per-condition Spearman cor(SAMac score, component score)
cors <- md %>% group_by(model, timepoint) %>%
  summarise(n = n(),
            cor_rev = cor(SAMac_sc, reversal_prog, method = "spearman"),
            cor_eme = cor(SAMac_sc, emergent_prog, method = "spearman"),
            p_rev = cor.test(SAMac_sc, reversal_prog, method = "spearman")$p.value,
            p_eme = cor.test(SAMac_sc, emergent_prog, method = "spearman")$p.value,
            .groups = "drop")
print(cors)
write.csv(cors, file.path(OUT, "mouse_samac_score_correlations.csv"), row.names = FALSE)

# (2) SAMac-high (top quartile within condition) vs rest: component scores
md2 <- md %>% group_by(model, timepoint) %>%
  mutate(samac_hi = SAMac_sc >= quantile(SAMac_sc, 0.75)) %>% ungroup()
hi_tests <- md2 %>% group_by(model, timepoint) %>%
  summarise(
    rev_hi = median(reversal_prog[samac_hi]), rev_lo = median(reversal_prog[!samac_hi]),
    eme_hi = median(emergent_prog[samac_hi]), eme_lo = median(emergent_prog[!samac_hi]),
    p_rev = wilcox.test(reversal_prog ~ samac_hi)$p.value,
    p_eme = wilcox.test(emergent_prog ~ samac_hi)$p.value, .groups = "drop")
print(hi_tests)
write.csv(hi_tests, file.path(OUT, "mouse_samac_hi_tests.csv"), row.names = FALSE)

# (3) does the SAMac-score--reversal correlation strengthen during regression? (Fisher z trend)
cors$z_rev <- atanh(cors$cor_rev); cors$z_eme <- atanh(cors$cor_eme)
print(cors %>% select(model, timepoint, cor_rev, cor_eme))
saveRDS(mo@meta.data, file.path(OUT, "mouse_metadata_with_samac_score.rds"))
cat("DONE\n")
