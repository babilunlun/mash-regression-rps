#!/usr/bin/env Rscript
# promote_fig2_gh.R — Fig2-S4 -> Fig2 主图 panel g；Fig2-S5 -> Fig2 主图 panel h
# 面板重出：去掉内部 a/b/c 小字母与 "Fig2-Sx" 总标题（主图统一标字母），内容/数据不变
# 然后重建 8 面板 Fig2_main
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
source("/workspace/fig1_theme.R")
suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(cowplot)})

T2  <- "/mnt/results/Figure2/tables"
F2  <- "/mnt/results/Figure2"
P2  <- "/workspace/refig/fig2_panels"
OUT <- "/workspace/refig/mains"
dir.create(P2, showWarnings = FALSE, recursive = TRUE)
comp_cols <- c(reversal = "#DCA7EB", emergent = "#BFDA8E", persistent = "#C8CACB", other = "#F0F0F0")
img <- function(p) wrap_elements(full = cowplot::ggdraw() + cowplot::draw_image(p))
tag_theme <- theme(plot.tag = element_text(size = 8, face = "bold", family = "Liberation Sans"))

## ============ panel g: GSE193066 独立验证（原 Fig2-S4 内容，无内部字母/总标题） ============
st2 <- read.csv(file.path(T2, "gse193066_staging_scores.csv"), stringsAsFactors = FALSE)
tr2 <- read.csv(file.path(T2, "gse193066_staging_trend.csv"), stringsAsFactors = FALSE)
dlt <- read.csv(file.path(T2, "gse193066_paired_deltas.csv"), stringsAsFactors = FALSE)
tst <- read.csv(file.path(T2, "gse193066_paired_tests.csv"), stringsAsFactors = FALSE)
gc2 <- read.csv(file.path(T2, "gse193066_genelevel_concordance.csv"), stringsAsFactors = FALSE)

st2l <- st2 %>% pivot_longer(c(reversal, emergent), names_to = "component", values_to = "score") %>%
  mutate(component = factor(component, levels = c("reversal", "emergent")),
         fibrosis = factor(fibrosis, levels = 0:4, labels = paste0("F", 0:4)))
tr2l <- tr2 %>% mutate(component = factor(geneset, levels = c("reversal", "emergent")),
                       label = sprintf("rho = %.2f\nFDR = %.3f", rho, fdr))
s4a <- ggplot(st2l, aes(fibrosis, score)) +
  geom_boxplot(aes(fill = component), alpha = 0.75, outlier.size = 0.7, outlier.alpha = 0.5,
               linewidth = 0.35, width = 0.62) +
  geom_line(stat = "summary", fun = median, aes(group = 1), color = "grey20", linewidth = 0.8) +
  geom_point(stat = "summary", fun = median, color = "grey20", size = 1.8) +
  geom_text(data = tr2l, aes(x = 4.6, y = -Inf, label = label), hjust = 1, vjust = -0.4,
            size = 3, lineheight = 1.05, inherit.aes = FALSE) +
  facet_wrap(~component, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = comp_cols, guide = "none") +
  labs(title = "Independent staging (1st biopsies, n=106)",
       x = "Fibrosis stage", y = "ssGSEA score") +
  theme_ipf_lf(11)

dlt2 <- dlt %>% pivot_longer(c(reversal, emergent), names_to = "component", values_to = "delta") %>%
  mutate(component = factor(component, levels = c("reversal", "emergent")),
         class = factor(class, levels = c("progressor", "stable", "regressor")))
tst2 <- tst %>% filter(geneset %in% c("reversal", "emergent")) %>%
  mutate(component = factor(geneset, levels = c("reversal", "emergent")),
         label = sprintf("regressor vs progressor\nWilcoxon P = %.3f", wilcox_p_reg_vs_prog))
class_cols2 <- c(progressor = "#E2EA4D", stable = "#C8CACB", regressor = "#B1C4D7")
s4b <- ggplot(dlt2, aes(class, delta, fill = class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, linewidth = 0.35, width = 0.55) +
  geom_point(position = position_jitter(width = 0.15, seed = 3), size = 1.1, alpha = 0.7) +
  geom_text(data = tst2, aes(x = 2, y = Inf, label = label), vjust = 1.3, size = 2.8,
            lineheight = 1.05, inherit.aes = FALSE) +
  facet_wrap(~component, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = class_cols2, name = NULL) +
  labs(title = "Paired deltas by histological course",
       subtitle = sprintf("gene-level concordance in regressors: %d/%d (%.0f%%), sign P = %.1e",
                          gc2$concordant, gc2$n, 100 * gc2$frac, gc2$sign_p),
       x = NULL, y = "Delta aligned score\n(+ = toward regression direction)") +
  theme_ipf_lf(11) + theme(legend.position = "none",
                           plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "#666666"))
p_g <- s4a + s4b + plot_layout(widths = c(1, 1.15))
save_panel(p_g, file.path(P2, "fig2g_gse193066"), width = 12.5, height = 4.4)

## ============ panel h: MASH 特异层（原 Fig2-S5 内容，无内部字母；GO 读已存表） ============
d <- read.csv(file.path(T2, "hcv_vs_mash_regression.csv"))
de <- d %>% filter(adj.P.Val_mash < 0.05) %>%
  mutate(concordant = sign(logFC_mash) == sign(logFC_hcv),
         cls = case_when(is_isg ~ "ISG (antiviral)",
                         concordant ~ "Shared (concordant in HCV)",
                         TRUE ~ "MASH-specific (discordant in HCV)"))
comp <- read.csv(file.path(T2, "components.csv")) %>% dplyr::select(gene, component)
de <- de %>% left_join(comp, by = "gene")
de$component[is.na(de$component)] <- "not in universe"
stopifnot(nrow(de) == 353)

set.seed(7)
cols <- c("Shared (concordant in HCV)" = "#C8CACB", "ISG (antiviral)" = "#B1C4D7",
          "MASH-specific (discordant in HCV)" = "#E2EA4D")
h_a <- ggplot(de, aes(x = logFC_mash, y = logFC_hcv)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80") +
  geom_point(data = de %>% filter(cls != "MASH-specific (discordant in HCV)"),
             aes(color = cls), size = 1.1, alpha = 0.75) +
  geom_point(data = de %>% filter(cls == "MASH-specific (discordant in HCV)"),
             aes(color = cls), size = 1.6, alpha = 0.9) +
  scale_color_manual(values = cols, name = NULL,
                     labels = c("ISG (antiviral)" = "ISG",
                                "MASH-specific (discordant in HCV)" = "MASH-specific",
                                "Shared (concordant in HCV)" = "Shared")) +
  theme_ipf_lf() + theme(legend.position = c(0.72, 0.85),
                         legend.text = element_text(size = 5.5),
                         legend.background = element_rect(fill = "white", color = NA),
                         plot.title = element_text(size = 7.5)) +
  labs(x = "MASH regression log2FC (GSE106737)", y = "HCV cure log2FC (GSE70779)") +
  ggtitle("Shared vs MASH-specific")

cmp <- de %>% count(component, cls) %>%
  group_by(cls) %>% mutate(frac = n / sum(n)) %>% ungroup()
cmp$component <- factor(cmp$component, levels = c("reversal", "emergent", "persistent", "other", "not in universe"))
n_p_ms <- sum(de$cls == "MASH-specific (discordant in HCV)" & de$component == "persistent")
n_ms   <- sum(de$cls == "MASH-specific (discordant in HCV)")
n_p_sh <- sum(de$cls == "Shared (concordant in HCV)" & de$component == "persistent")
n_sh   <- sum(de$cls == "Shared (concordant in HCV)")
ft <- fisher.test(matrix(c(n_p_ms, n_ms - n_p_ms, n_p_sh, n_sh - n_p_sh), nrow = 2))
cat(sprintf("persistent: MASH-specific %d/%d vs shared %d/%d, Fisher P=%.2e\n", n_p_ms, n_ms, n_p_sh, n_sh, ft$p.value))
stopifnot(n_ms == 75, n_sh == 255)
cmp$cls <- factor(cmp$cls, levels = c("MASH-specific (discordant in HCV)", "Shared (concordant in HCV)", "ISG (antiviral)"))
h_b <- ggplot(cmp, aes(x = cls, y = n, fill = component)) +
  geom_col(position = "fill", width = 0.65, alpha = 0.85) +
  scale_fill_manual(values = c(reversal = "#DCA7EB", emergent = "#BFDA8E", persistent = "#B1C4D7",
                               other = "#ECE9E2", `not in universe` = "#F5F5F5"),
                    name = "Component") +
  theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 6), legend.position = "bottom") +
  labs(x = NULL, y = "Fraction of genes") +
  annotate("text", x = 1, y = 1.06, label = sprintf("persistent enrichment vs shared: Fisher P=%.1e", ft$p.value),
           size = 2.2, hjust = 0.35, family = "Liberation Sans") +
  coord_cartesian(clip = "off") +
  ggtitle("Component membership")

eg <- read.csv(file.path(T2, "mash_specific_go_bp.csv"), stringsAsFactors = FALSE)
min_fdr <- min(eg$p.adjust)
eg <- eg %>% arrange(pvalue) %>% head(10)
eg$Description <- factor(eg$Description, levels = rev(eg$Description))
fdr_note <- if (min_fdr >= 0.05) " (nominal; no term survives FDR)" else ""
h_c <- ggplot(eg, aes(x = -log10(pvalue), y = Description, fill = Count)) +
  geom_col(alpha = 0.85) +
  scale_fill_gradient(low = "#B1C4D7", high = "#5B7FA6") +
  theme_ipf_lf() + theme(legend.position = "bottom",
                         axis.text.y = element_text(size = 5.8),
                         plot.title = element_text(size = 7.5)) +
  labs(x = expression(-log[10](italic(P)[nominal])), y = NULL) +
  ggtitle(paste0("GO-BP, MASH-specific", fdr_note))

p_h <- h_a + h_b + h_c + plot_layout(widths = c(1.15, 0.75, 1.1))
save_panel(p_h, file.path(P2, "fig2h_mash_specific"), width = 11.5, height = 3.8)

## ============ Fig2 主图重建：8 面板（a RRHO | b 分量 | c 分期 | d 跨物种 | e 小鼠时程 | f HCV | g GSE193066 | h MASH特异层） ============
# a–f 单元格保持原 6 面板版尺寸(3.6x2.6in)；g/h 通栏行高按面板原始宽高比精确匹配
fig2 <- img(file.path(F2, "Fig2a_RRHO.png")) + img(file.path(F2, "Fig2b_components.png")) +
        img(file.path(P2, "fig2c_staging.png")) + img(file.path(F2, "Fig2d_cross_species.png")) +
        img(file.path(P2, "fig2e_mouse_temporal.png")) + img(file.path(F2, "Fig2f_HCV_concordance.png")) +
        img(file.path(P2, "fig2g_gse193066.png")) + img(file.path(P2, "fig2h_mash_specific.png")) +
  plot_layout(design = "AAABBB\nAAABBB\nCCCDDD\nCCCDDD\nEEEFFF\nEEEFFF\nGGGGGG\nHHHHHH",
              heights = c(1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 2.53, 2.38)) +
  plot_annotation(tag_levels = "a") & tag_theme
ggsave(file.path(OUT, "Fig2_main.svg"), fig2, width = 7.2, height = 12.71, device = svglite::svglite)
ggsave(file.path(OUT, "Fig2_main.png"), fig2, width = 7.2, height = 12.71, dpi = 600)
cat("FIG2_PROMOTE_DONE\n")
