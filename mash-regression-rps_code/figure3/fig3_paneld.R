#!/usr/bin/env Rscript
# Figure 3 panel d: CellChat regression-specific interactions
suppressMessages({ library(dplyr); library(ggplot2); library(patchwork); library(tidyr) })
source("/workspace/fig1_theme.R")
OUT <- "/workspace/fig3_out"; FIG <- "/workspace/fig3_fig"

## (left) pathway deltas: shared up/down in both models
pd <- read.csv(file.path(OUT, "cellchat_pathway_deltas.csv"))
# normalize by total pathway strength per condition (remove global scale differences)
pd <- pd %>% group_by(model) %>%
  mutate(rel_peak = peak / sum(peak), rel_reg = reg6wk / sum(reg6wk),
         delta = rel_reg - rel_peak) %>% ungroup()
pw <- pd %>% select(model, pathway, delta) %>%
  pivot_wider(names_from = model, values_from = delta, values_fill = 0) %>%
  rename(delta_CCl4 = CCl4, `delta_FAT-MASH` = `FAT-MASH`)
pw <- pw %>% mutate(dir_ccl4 = sign(delta_CCl4), dir_fat = sign(`delta_FAT-MASH`),
                    shared = ifelse(dir_ccl4 == dir_fat, ifelse(dir_ccl4 > 0, "up_both", "down_both"), "divergent"),
                    mean_delta = (delta_CCl4 + `delta_FAT-MASH`) / 2)
top <- pw %>% filter(shared != "divergent") %>%
  group_by(shared) %>% slice_max(abs(mean_delta), n = 8, with_ties = FALSE) %>% ungroup()
top$pathway <- factor(top$pathway, levels = top$pathway[order(top$mean_delta)])
pd_l <- pw %>% filter(pathway %in% top$pathway) %>%
  select(pathway, delta_CCl4, `delta_FAT-MASH`) %>%
  pivot_longer(-pathway, names_to = "model", values_to = "delta")
pd_l$model <- sub("delta_", "", pd_l$model)
pd_l$pathway <- factor(pd_l$pathway, levels = levels(top$pathway))
pd1 <- ggplot(pd_l, aes(x = pathway, y = delta, fill = model)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  coord_flip() + geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(CCl4 = palette8[2], `FAT-MASH` = palette8[4])) +
  theme_ipf_lf() + labs(x = NULL, y = "Relative pathway strength delta (reg6wk - peak)", fill = NULL)

## (right) myeloid->HSC gained LR pairs in regression
mye_gain <- bind_rows(lapply(c("CCl4", "FAT-MASH"), function(mdl) {
  pk <- read.csv(file.path(OUT, paste0("cellchat_myeloid_to_hsc_", mdl, "_peak.csv")))
  rg <- read.csv(file.path(OUT, paste0("cellchat_myeloid_to_hsc_", mdl, "_reg6wk.csv")))
  key <- function(d) paste(d$source, d$ligand, d$receptor, sep = "|")
  gained <- rg[!key(rg) %in% key(pk), ]
  gained$model <- mdl
  gained
}))
# keep pairs present with prob >= 0.005, top by prob
mg <- mye_gain %>% filter(prob >= 0.005) %>%
  mutate(lr = paste0(ligand, "-", receptor))
# order by mean prob
ord <- mg %>% group_by(lr) %>% summarise(m = mean(prob), .groups = "drop") %>% arrange(desc(m))
mg$lr <- factor(mg$lr, levels = rev(ord$lr))
pd2 <- ggplot(mg, aes(x = model, y = lr, size = prob, color = source)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c(`SAMac (Trem2+)` = "#C7B0C7", `Monocyte/MF` = "#A7C0DF",
                                Kupffer = "#BEF3C5", Neutrophil = "#FDB863", DC = "#F4A582")) +
  theme_ipf_lf() + theme(axis.text.y = element_text(size = 6)) +
  labs(x = NULL, y = NULL, size = "Comm. prob.", color = "Source") +
  ggtitle("Myeloid-to-HSC interactions gained in regression")
pd <- pd1 | pd2
save_panel(pd, file.path(FIG, "fig3d_cellchat"), width = 12, height = 5.5)
cat("PANEL D DONE\n")
