#!/usr/bin/env Rscript
# Figure 3 panel f: MISTy neighborhood importances (mouse Visium) + CosMx LR colocalization (human)
suppressMessages({ library(dplyr); library(ggplot2); library(patchwork); library(tidyr) })
source("/workspace/fig1_theme.R")
OUT <- "/workspace/fig3_out"; FIG <- "/workspace/fig3_fig"

## (left) MISTy: top neighborhood (paraview) predictors of reversal score, per slide
imp <- bind_rows(lapply(c("peak", "regression"), function(sn) {
  read.csv(file.path(OUT, paste0("misty_importances_", sn, ".csv"))) %>% mutate(slide = sn)
}))
mi <- imp %>% filter(Target == "reversal_prog", view == "para.200",
                     !Predictor %in% c("reversal_prog", "emergent_prog")) %>%
  group_by(slide) %>% slice_max(Importance, n = 8, with_ties = FALSE) %>% ungroup()
mi$Predictor <- dplyr::recode(mi$Predictor,
  `Endothelial.LSEC` = "Endothelial/LSEC", `T.NK` = "T/NK", `B.cell` = "B cell",
  `Monocyte.MF` = "Monocyte/MF", `HSC.myofibroblast` = "HSC/myofibroblast",
  `SAMac..Trem2..` = "SAMac (Trem2+)", .default = mi$Predictor)
ord <- mi %>% group_by(Predictor) %>% summarise(m = mean(Importance), .groups = "drop") %>% arrange(m)
mi$Predictor <- factor(mi$Predictor, levels = ord$Predictor)
pf1 <- ggplot(mi, aes(x = Predictor, y = Importance, fill = slide)) +
  geom_col(position = position_dodge(0.7), width = 0.65) + coord_flip() +
  scale_fill_manual(values = c(peak = "#C8CACB", regression = "#BFDA8E")) +
  theme_ipf_lf() + labs(x = NULL, y = "MISTy neighborhood importance", fill = NULL) +
  ggtitle("Neighborhood predictors of reversal score (Visium)")

## (right) CosMx LR colocalization
lr <- read.csv(file.path(OUT, "cosmx_lr_colocalization.csv"))
lr$sig <- !is.na(lr$fdr) & lr$fdr < 0.05
lr$lr <- factor(lr$lr, levels = lr %>% group_by(lr) %>% summarise(m = mean(log2fc), .groups = "drop") %>%
                  arrange(m) %>% pull(lr))
# manual dodge to keep stars aligned after coord_flip
lr$ypos <- as.numeric(lr$lr) + ifelse(lr$slide == "Slide1", -0.18, 0.18)
pf2 <- ggplot(lr, aes(y = ypos, x = log2fc, color = slide)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(size = 2) +
  geom_text(data = lr %>% filter(sig), aes(x = log2fc + 0.06 * sign(log2fc)), label = "*",
            size = 5, show.legend = FALSE) +
  scale_y_continuous(breaks = seq_along(levels(lr$lr)), labels = levels(lr$lr)) +
  scale_color_manual(values = c(Slide1 = "#DCA7EB", Slide2 = "#A7C0DF")) +
  theme_ipf_lf() + labs(y = NULL, x = "Receptor log2FC (near vs far stellate)", color = NULL) +
  ggtitle("Ligand+ macrophage proximity effect (CosMx)")
pf <- pf1 | pf2
save_panel(pf, file.path(FIG, "fig3f_misty_cosmx"), width = 11, height = 4.5)
cat("PANEL F DONE\n")
