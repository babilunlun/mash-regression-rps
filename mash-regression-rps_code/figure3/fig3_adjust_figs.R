#!/usr/bin/env Rscript
# Phase-1 adjustments: fig3S6 (SAMac continuous analysis) + fig3b v2 (FDR stars)
suppressMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(tidyr) })
source("/workspace/fig1_theme.R")
OUT <- "/workspace/fig3_out"; FIG <- "/workspace/fig3_fig"

## ---------- fig3S6: SAMac continuous ----------
cors <- read.csv(file.path(OUT, "mouse_samac_score_correlations.csv"))
cors$timepoint <- factor(cors$timepoint, levels = c("peak", "reg2wk", "reg6wk"))
cl <- cors %>% select(model, timepoint, cor_rev, cor_eme) %>%
  pivot_longer(c(cor_rev, cor_eme), names_to = "comp", values_to = "rho")
cl$comp <- sub("cor_", "", cl$comp)
s6a <- ggplot(cl, aes(x = timepoint, y = rho, fill = comp)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  facet_wrap(~model) +
  scale_fill_manual(values = c(rev = "#DCA7EB", eme = "#BFDA8E"),
                    breaks = c("rev", "eme"), labels = c("reversal", "emergent")) +
  theme_ipf_lf() + ylim(0, 0.75) +
  labs(x = NULL, y = "Spearman rho (SAMac score vs component)", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

hi <- read.csv(file.path(OUT, "mouse_samac_hi_tests.csv"))
hi$timepoint <- factor(hi$timepoint, levels = c("peak", "reg2wk", "reg6wk"))
hi_l <- hi %>% select(model, timepoint, rev_hi, rev_lo) %>%
  pivot_longer(c(rev_hi, rev_lo), names_to = "group", values_to = "med")
hi_l$group <- factor(hi_l$group, levels = c("rev_lo", "rev_hi"), labels = c("SAMac-low", "SAMac-high"))
s6b <- ggplot(hi_l, aes(x = timepoint, y = med, fill = group)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  facet_wrap(~model) +
  scale_fill_manual(values = c(`SAMac-low` = "#C8CACB", `SAMac-high` = "#C7B0C7")) +
  theme_ipf_lf() + labs(x = NULL, y = "Median reversal UCell", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
s6 <- s6a | s6b
save_panel(s6, file.path(FIG, "fig3S6_samac_continuous"), width = 9, height = 3.8)

## ---------- fig3b v2: add FDR stars ----------
h <- readRDS(file.path(OUT, "human_sc_scored.rds"))
mo <- readRDS(file.path(OUT, "mouse_sc_scored.rds"))
mk_long <- function(obj, group_var, species) {
  md <- obj@meta.data
  bind_rows(lapply(c("reversal", "emergent"), function(comp) {
    data.frame(species = species, comp = comp, celltype = md$celltype_label,
               grp = md[[group_var]], score = md[[paste0(comp, "_prog")]])
  }))
}
dl <- bind_rows(mk_long(h, "group", "Human"), mk_long(mo, "grp2", "Mouse"))
dl <- dl %>% filter(!celltype %in% c("Other", "Erythrocyte", "Platelet"))
dl$celltype <- factor(dl$celltype)
dl$grp <- factor(dl$grp, levels = c("Healthy", "MASLD", "peak", "reg2wk", "reg6wk"))
grp_pal <- c(Healthy = "#C8CACB", MASLD = "#E2EA4D", peak = "#E2EA4D", reg2wk = "#B1C4D7", reg6wk = "#BFDA8E")

# stars from test tables
th <- read.csv(file.path(OUT, "ucell_tests_human.csv")) %>%
  filter(!is.na(fdr)) %>% transmute(species = "Human", comp, celltype, fdr,
                                    grp = "MASLD", y = pmax(med_healthy, med_masld))
tm <- read.csv(file.path(OUT, "ucell_tests_mouse.csv")) %>%
  filter(!is.na(fdr)) %>% transmute(species = "Mouse", comp, celltype, fdr,
                                    grp = "reg6wk", y = pmax(med_peak, med_reg6wk))
stars <- bind_rows(th, tm) %>%
  mutate(star = case_when(fdr < 0.001 ~ "***", fdr < 0.01 ~ "**", fdr < 0.05 ~ "*", TRUE ~ "")) %>%
  filter(star != "") %>%
  filter(celltype %in% unique(dl$celltype))
# y position: per facet max score * 1.05
ymax <- dl %>% filter(grp %in% c("Healthy", "MASLD", "peak", "reg6wk")) %>%
  group_by(species, comp, celltype) %>% summarise(ymax = max(score), .groups = "drop")
stars <- stars %>% left_join(ymax, by = c("species", "comp", "celltype")) %>%
  mutate(ypos = ymax * 1.06)

pb <- ggplot(dl %>% filter(grp %in% c("Healthy", "MASLD", "peak", "reg6wk")),
             aes(x = celltype, y = score, fill = grp)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.2, linewidth = 0.2, position = position_dodge(0.75)) +
  facet_grid(comp ~ species, scales = "free") +
  scale_fill_manual(values = grp_pal) +
  geom_text(data = stars, aes(x = celltype, y = ypos, label = star, group = grp),
            inherit.aes = FALSE, size = 2.6, vjust = 0.3) +
  theme_ipf_lf() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  labs(x = NULL, y = "UCell score", fill = NULL)
save_panel(pb, file.path(FIG, "fig3b_ucell_scores"), width = 11, height = 6)
cat("ADJUST FIGS DONE\n")
