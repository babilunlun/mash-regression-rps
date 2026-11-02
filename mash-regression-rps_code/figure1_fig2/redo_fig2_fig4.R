# redo_fig2_fig4.R — 仅布局/字母修正
# Fig2: 重出 staging(新c) 与 mouse temporal(新e) 两块面板 —— 设计与原版完全一致，仅去掉标题中烤入的旧字母前缀
# Fig4: 改为通栏堆叠布局（a Visium / b GeoMx / c MISTy | d CosMx），修正 MISTy 标题裁切
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
source("/workspace/fig1_theme.R")
suppressMessages({library(data.table); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(cowplot)})

OUT <- "/workspace/refig/mains"
P2  <- "/workspace/refig/fig2_panels"
dir.create(P2, showWarnings = FALSE, recursive = TRUE)

img <- function(p) wrap_elements(full = cowplot::ggdraw() + cowplot::draw_image(p))
tag_theme <- theme(plot.tag = element_text(size = 8, face = "bold", family = "Liberation Sans"))
save_comb <- function(plot, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".svg")), plot, width = w, height = h, device = svglite::svglite)
  ggsave(file.path(OUT, paste0(name, ".png")), plot, width = w, height = h, dpi = 600)
}

T2 <- "/mnt/results/Figure2/tables"
F2 <- "/mnt/results/Figure2"
F3 <- "/mnt/results/Figure3"
comp_cols <- c(reversal = "#DCA7EB", emergent = "#BFDA8E", persistent = "#C8CACB", other = "#F0F0F0")

## ============ Fig2 新c：human staging（原 Fig2e，同设计，标题去字母） ============
st_scores <- read.csv(file.path(T2, "staging_scores_GSE135251.csv"), stringsAsFactors = FALSE)
st_trend  <- read.csv(file.path(T2, "staging_trend_GSE135251.csv"), stringsAsFactors = FALSE)
st_long <- st_scores %>%
  pivot_longer(c(reversal, emergent), names_to = "component", values_to = "score") %>%
  mutate(component = factor(component, levels = c("reversal", "emergent")),
         fibrosis = factor(fibrosis, levels = 0:4, labels = paste0("F", 0:4)))
trend_lab <- st_trend %>% filter(geneset %in% c("reversal", "emergent")) %>%
  mutate(component = factor(geneset, levels = c("reversal", "emergent")),
         label = sprintf("Spearman rho = %.2f\nFDR = %.1e", rho, fdr))

p_staging <- ggplot(st_long, aes(fibrosis, score)) +
  geom_boxplot(aes(fill = component), alpha = 0.75, outlier.size = 0.7, outlier.alpha = 0.5,
               linewidth = 0.35, width = 0.62) +
  geom_line(stat = "summary", fun = median, aes(group = 1), color = "grey20", linewidth = 0.8) +
  geom_point(stat = "summary", fun = median, color = "grey20", size = 1.8) +
  geom_text(data = trend_lab, aes(x = 4.6, y = -Inf, label = label), hjust = 1, vjust = -0.4,
            size = 3, lineheight = 1.05, inherit.aes = FALSE) +
  facet_wrap(~component, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = comp_cols, guide = "none") +
  labs(title = "Component scores across human fibrosis stages (GSE135251)",
       x = "Fibrosis stage", y = "ssGSEA score") +
  theme_ipf_lf(11)
save_panel(p_staging, file.path(P2, "fig2c_staging"), width = 8, height = 4.2)

## ============ Fig2 新e：mouse temporal（原 Fig2c，同设计，标题去字母） ============
traj_mean <- read.csv(file.path(T2, "trajectory_means_GSE294630.csv"), stringsAsFactors = FALSE)
comp_orth <- read.csv(file.path(T2, "component_orthologs.csv"), stringsAsFactors = FALSE)
el_sema   <- read.csv(file.path(T2, "component_earlylate_sema.csv"), stringsAsFactors = FALSE)
el_chow   <- read.csv(file.path(T2, "component_earlylate_chowsw.csv"), stringsAsFactors = FALSE)

grps <- c("baseline", "veh12", "veh18", "chowsw12", "chowsw18", "sema12", "sema18")
tm <- traj_mean %>% filter(ensg %in% comp_orth$ensembl)
tm <- tm[match(comp_orth$ensembl, tm$ensg), ]
stopifnot(all(tm$ensg == comp_orth$ensembl))
z <- as.matrix(tm[, grps]); z <- t(scale(t(z)))
sgn <- ifelse(comp_orth$regr_logFC > 0, 1, -1)
za <- sweep(z, 1, sgn, `*`)
weeks <- c(baseline = 6, veh12 = 12, veh18 = 18, chowsw12 = 12, chowsw18 = 18, sema12 = 12, sema18 = 18)
arms <- list(vehicle = c("baseline", "veh12", "veh18"),
             `chow-switch` = c("baseline", "chowsw12", "chowsw18"),
             semaglutide = c("baseline", "sema12", "sema18"))
traj_df <- bind_rows(lapply(names(arms), function(arm) {
  bind_rows(lapply(c("reversal", "emergent"), function(cp) {
    idx <- which(comp_orth$component == cp)
    m <- colMeans(za[idx, arms[[arm]], drop = FALSE])
    se <- apply(za[idx, arms[[arm]], drop = FALSE], 2, function(v) sd(v) / sqrt(length(v)))
    data.frame(component = cp, arm = arm, week = weeks[arms[[arm]]], mean = m, se = se)
  }))
}))
traj_df$component <- factor(traj_df$component, levels = c("reversal", "emergent"))
arm_cols <- c(vehicle = "#C8CACB", `chow-switch` = "#B1C4D7", semaglutide = "#E2EA4D")

c1 <- ggplot(traj_df, aes(week, mean, color = arm, group = arm)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.5, linewidth = 0.4) +
  facet_wrap(~component, nrow = 1) +
  scale_color_manual(values = arm_cols, name = NULL) +
  scale_x_continuous(breaks = c(6, 12, 18)) +
  labs(title = "Component trajectories in mouse regression (GSE294630)",
       x = "Weeks on CDAA-HFD (treatment from week 6)",
       y = "Direction-aligned z-score\n(+ = regression-favored)") +
  theme_ipf_lf(11) +
  theme(legend.position = c(0.03, 0.97), legend.justification = c(0, 1),
        legend.text = element_text(size = 8.5),
        legend.background = element_rect(fill = alpha("white", 0.75), color = NA))

el <- bind_rows(el_sema, el_chow) %>%
  mutate(component = factor(component, levels = c("reversal", "emergent")),
         class = factor(class, levels = c("sustained", "early_only", "late_only", "none")),
         arm = ifelse(arm == "sema", "semaglutide", "chow-switch")) %>%
  group_by(arm, component) %>% mutate(frac = n / sum(n)) %>% ungroup()
class_cols <- c(sustained = "#A7C0DF", early_only = "#BEF3C5", late_only = "#E2EA4D", none = "#EFEFEF")

c2 <- ggplot(el, aes(component, frac, fill = class)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(data = el %>% filter(frac > 0.06),
            aes(label = n), position = position_stack(vjust = 0.5), size = 2.6) +
  facet_wrap(~arm, nrow = 1) +
  scale_fill_manual(values = class_cols, name = NULL,
                    labels = c("sustained (12w+18w)", "early only (12w)", "late only (18w)", "not DE")) +
  labs(title = "   Temporal class of component orthologs", x = NULL, y = "Fraction of genes") +
  theme_ipf_lf(11) + theme(legend.position = "right", legend.text = element_text(size = 8.5))

p_temporal <- c1 + c2 + plot_layout(widths = c(1.25, 1))
save_panel(p_temporal, file.path(P2, "fig2e_mouse_temporal"), width = 11, height = 4.2)

## ============ Fig 2 重新拼装（a RRHO | b 分量 | c 分期 | d 跨物种 | e 小鼠时程 | f HCV） ============
fig2 <- img(file.path(F2, "Fig2a_RRHO.png")) + img(file.path(F2, "Fig2b_components.png")) +
        img(file.path(P2, "fig2c_staging.png")) + img(file.path(F2, "Fig2d_cross_species.png")) +
        img(file.path(P2, "fig2e_mouse_temporal.png")) + img(file.path(F2, "Fig2f_HCV_concordance.png")) +
  plot_layout(design = "AAABBB\nAAABBB\nCCCDDD\nCCCDDD\nEEEFFF\nEEEFFF") +
  plot_annotation(tag_levels = "a") & tag_theme
save_comb(fig2, "Fig2_main", 7.2, 7.8)

## ============ Fig 4：MISTy/CosMx 标题定位修正 + 通栏堆叠布局 ============
lab_map <- c(Endothelial_LSEC = "Endothelial/LSEC", T_NK = "T/NK", pDC = "pDC", B.cell = "B cell",
             Kupffer = "Kupffer", Neutrophil = "Neutrophil", Monocyte_MF = "Monocyte/MF",
             Erythrocyte = "Erythrocyte", `SAMac..Trem2..` = "SAMac (Trem2+)",
             HSC_myofibroblast = "HSC/myofibroblast")
get_imp <- function(f) {
  d <- fread(f)
  s <- d[view == "para.200" & Target == "reversal_prog" & Predictor %in% names(lab_map),
         .(Predictor, Importance)]
  s[, ct := lab_map[Predictor]]
  s[, .(ct, Importance)]
}
pk <- get_imp(file.path(F3, "tables/misty_importances_peak.csv"))
rg <- get_imp(file.path(F3, "tables/misty_importances_regression.csv"))
mm <- merge(pk[, .(ct, peak = Importance)], rg[, .(ct, regression = Importance)], by = "ct")
mm[, ct := factor(ct, levels = mm$ct[order(mm$regression)])]
ml <- melt(mm, id.vars = "ct", variable.name = "tp", value.name = "imp")

pMISTy <- ggplot(ml, aes(x = imp, y = ct)) +
  geom_col(data = ml[tp == "peak"], fill = "#C8CACB", alpha = 0.85, width = 0.68) +
  geom_col(data = ml[tp == "regression"], fill = "#BFDA8E", alpha = 0.85, width = 0.68) +
  labs(x = "MISTy neighborhood importance", y = NULL,
       title = "Neighborhood predictors of reversal score (Visium)") +
  theme_ipf_lf(base = 9) +
  theme(legend.position = "none",
        plot.title = element_text(size = 9.5), plot.title.position = "plot")

cx <- fread(file.path(F3, "tables/cosmx_lr_colocalization.csv"))
lr_ord <- cx[, .(m = mean(log2fc)), by = lr][order(-m)]$lr
cx[, lr := factor(lr, levels = rev(lr_ord))]
pCosMx <- ggplot(cx, aes(x = log2fc, y = lr, color = slide)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_text(data = cx[fdr < 0.05], aes(label = "*"), size = 5, nudge_x = 0.045,
            vjust = 0.3, show.legend = FALSE, family = "Liberation Sans") +
  scale_color_manual(values = c(Slide1 = "#DCA7EB", Slide2 = "#A7C0DF"), name = NULL) +
  labs(x = "Receptor log2FC (near vs far stellate)", y = NULL,
       title = "Ligand+ macrophage proximity effect (CosMx)") +
  theme_ipf_lf(base = 9) +
  theme(plot.title = element_text(size = 9.5), plot.title.position = "plot",
        legend.position = c(0.88, 0.85))

save_panel(pMISTy, file.path(OUT, "fig4c_misty"), width = 4.6, height = 3.4)
save_panel(pCosMx, file.path(OUT, "fig4d_cosmx"), width = 4.6, height = 3.4)

fig4 <- img(file.path(F3, "fig3e_mouse_visium.png")) + img(file.path(F3, "fig3e_human_geomx.png")) +
        wrap_elements(full = pMISTy) + wrap_elements(full = pCosMx) +
  plot_layout(design = "AAAAAAAA\nAAAAAAAA\nAAAAAAAA\nBBBBBBBB\nBBBBBBBB\nBBBBBBBB\nCCCCDDDD\nCCCCDDDD\nCCCCDDDD") +
  plot_annotation(tag_levels = "a") & tag_theme
save_comb(fig4, "Fig4_main", 7.2, 9.0)

cat("REDO_DONE\n")
