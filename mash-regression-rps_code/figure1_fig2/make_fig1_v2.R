# make_fig1_v2.R — Fig 1 设计重做（5 版：a banner + b 火山 + c GO 双向点图 + d 配对斜率 + e 双层散点）
# 数据：/mnt/results/Figure1/tables（归档表，统计不变，纯设计层重绘）
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
source("/workspace/fig1_theme.R")

out <- "/workspace/refig/fig1"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
tab <- "/mnt/results/Figure1/tables"

de    <- fread(file.path(tab, "de_paired41.csv"))
mod   <- fread(file.path(tab, "haas2019_modules.csv"))
rep   <- fread(file.path(tab, "replication_concordance.csv"))
go_dn <- fread(file.path(tab, "enrich_GOBP_paired41_down.csv"))
go_up <- fread(file.path(tab, "enrich_GOBP_paired41_up.csv"))
ssg   <- fread(file.path(tab, "ssgsea_blue_GSE106737.csv"))

blue <- mod[Module == "blue", GeneSymbols]

# 颜色语义（全图一致）：下调/免疫回撤 = #A7C0DF；上调/代谢恢复 = #BFDA8E；中性灰 = #C8CACB；警示 = #DCA7EB
C_DN <- "#A7C0DF"; C_UP <- "#BFDA8E"; C_NS <- "#C8CACB"; C_WARN <- "#DCA7EB"

## ---------------- Panel b：火山图（Haas blue 叠加 + 直接标注） ----------------
de[, `:=`(is_blue = gene %in% blue, sig = adj.P.Val < 0.05)]
de[, grp := fcase(
  !sig, "ns",
  sig & logFC < 0 & is_blue, "down_blue",
  sig & logFC < 0, "down",
  sig & logFC > 0 & is_blue, "up_blue",
  sig & logFC > 0, "up")]
print(de[, .N, by = grp])

# Fisher：blue 模块在下调 DE 中的富集（核对 README: 107/250, OR=32.3）
n_univ <- nrow(de); n_blue <- sum(de$is_blue)
n_dn <- nrow(de[grp %in% c("down", "down_blue")]); n_ov <- nrow(de[grp == "down_blue"])
ft <- fisher.test(matrix(c(n_ov, n_dn - n_ov, n_blue - n_ov, n_univ - n_blue - (n_dn - n_ov)), 2, 2))
cat(sprintf("blue overlap %d/%d, OR=%.1f, p=%.2e\n", n_ov, n_dn, ft$estimate, ft$p.value))

thr_y <- -log10(max(de$P.Value[de$sig]))
lab_dn <- de[grp %in% c("down", "down_blue")][order(P.Value)][1:8]
lab_up <- de[grp %in% c("up", "up_blue")][order(P.Value)][1:4]

pB <- ggplot(de, aes(logFC, -log10(P.Value))) +
  geom_point(data = de[grp == "ns"], color = C_NS, alpha = 0.25, size = 0.5) +
  geom_point(data = de[grp == "down"], color = C_DN, alpha = 0.45, size = 0.9) +
  geom_point(data = de[grp == "up"], color = C_UP, alpha = 0.75, size = 0.9) +
  geom_point(data = de[grp == "down_blue"],
             shape = 21, fill = C_DN, color = "grey25", stroke = 0.25, size = 1.5) +
  geom_point(data = de[grp == "up_blue"],
             shape = 21, fill = C_UP, color = "grey25", stroke = 0.25, size = 1.5) +
  geom_hline(yintercept = thr_y, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_text_repel(data = rbind(lab_dn, lab_up), aes(label = gene),
                  size = 2.2, family = "Liberation Sans", fontface = "italic",
                  segment.size = 0.25, segment.color = "grey50",
                  box.padding = 0.35, point.padding = 0.2, max.overlaps = 30,
                  min.segment.length = 0.2, seed = 42) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.05, vjust = -0.7, size = 2.4,
           label = sprintf("%d down (FDR<0.05)", n_dn), color = "grey20") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.7, size = 2.4,
           label = sprintf("%d up (FDR<0.05)", nrow(de[grp %in% c("up", "up_blue")])), color = "grey20") +
  annotate("text", x = -Inf, y = thr_y, hjust = -0.05, vjust = 1.3, size = 2.2, color = "grey30",
           label = sprintf("%d/%d down-DE in Haas blue module\nOR = %.1f, P < 2e-16", n_ov, n_dn, ft$estimate)) +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.12))) +
  labs(x = "log2 fold change (1-year vs baseline)", y = "-log10 P") +
  theme_ipf_lf(base = 7) +
  theme(plot.margin = margin(5, 8, 5, 5))

## ---------------- Panel c：GO 双向点图 ----------------
n_dn_term <- 4; n_up_term <- 3
g <- rbind(
  head(go_dn[order(p.adjust)], n_dn_term)[, .(Description, Count, p.adjust, dir = "down")],
  head(go_up[order(p.adjust)], n_up_term)[, .(Description, Count, p.adjust, dir = "up")])
g[, x := -log10(p.adjust) * fifelse(dir == "down", -1, 1)]
g[, term := sapply(Description, function(s) paste(strwrap(s, 32), collapse = "\n"))]
g[, term := factor(term, levels = g$term[order(g$x)])]
ytop <- length(levels(g$term)) + 0.6

pC <- ggplot(g, aes(x = x, y = term)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = x, yend = term), color = "grey75", linewidth = 0.4) +
  geom_point(aes(size = Count, fill = dir), shape = 21, color = "grey30", stroke = 0.25, alpha = 0.9) +
  scale_fill_manual(values = c(down = C_DN, up = C_UP), guide = "none") +
  scale_size_continuous(range = c(1.8, 4.5), guide = "none") +
  scale_y_discrete(expand = expansion(add = c(0.8, 0.8))) +
  labs(x = "\u2190 down in regression      signed -log10 FDR      up in regression \u2192", y = NULL) +
  theme_ipf_lf(base = 7) +
  theme(axis.text.y = element_text(size = 5.2, lineheight = 0.82),
        plot.margin = margin(5, 8, 5, 5))

## ---------------- Panel d：模块配对斜率图 ----------------
pr <- ssg[sample_type == "paired"]
pr[, timepoint := factor(timepoint, levels = c("baseline", "followup"))]
pw <- dcast(pr, pair_id ~ timepoint, value.var = "blue_score")
pw[, delta := followup - baseline]
n_dec <- sum(pw$delta < 0); med_d <- median(pw$delta)
cat(sprintf("paired: %d/%d decreased, median delta = %.3f\n", n_dec, nrow(pw), med_d))
pr <- merge(pr, pw[, .(pair_id, delta)], by = "pair_id")
med <- pr[, .(blue_score = median(blue_score)), by = timepoint]

pD <- ggplot(pr, aes(timepoint, blue_score)) +
  geom_line(aes(group = pair_id, color = delta < 0), alpha = 0.55, linewidth = 0.4) +
  geom_point(aes(color = delta < 0), size = 0.9, alpha = 0.7) +
  geom_line(data = med, aes(group = 1), color = "black", linewidth = 1.1) +
  geom_point(data = med, color = "black", size = 2) +
  scale_color_manual(values = c(`TRUE` = C_DN, `FALSE` = C_WARN), guide = "none") +
  annotate("text", x = 1.5, y = Inf, vjust = 1.5, size = 2.4,
           label = sprintf("%d/%d decreased\nmedian delta = %.2f", n_dec, nrow(pw), med_d)) +
  scale_x_discrete(labels = c("Baseline", "1-year\nfollow-up")) +
  labs(x = NULL, y = "NASH module score (ssGSEA)") +
  theme_ipf_lf(base = 7) +
  theme(plot.margin = margin(5, 8, 5, 5))

## ---------------- Panel e：独立重复双层散点 ----------------
rep[, de1 := adj.P.Val_g106737 < 0.05]
rep[, conc := sign(logFC_g106737) == sign(logFC_g83452)]
de_genes <- rep[de1 == TRUE]
n_conc <- sum(de_genes$conc); n_de <- nrow(de_genes)
cat(sprintf("replication: %d/%d concordant (%.1f%%)\n", n_conc, n_de, 100 * n_conc / n_de))

pE <- ggplot(rep, aes(logFC_g106737, logFC_g83452)) +
  geom_point(data = rep[de1 == FALSE], color = C_NS, alpha = 0.12, size = 0.4) +
  geom_point(data = de_genes[conc == TRUE], aes(color = logFC_g106737 < 0), size = 1.1, alpha = 0.85) +
  geom_point(data = de_genes[conc == FALSE], color = C_WARN, size = 2, shape = 21, stroke = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  scale_color_manual(values = c(`TRUE` = C_DN, `FALSE` = C_UP), guide = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.4, size = 2.2,
           label = sprintf("Spearman \u03c1 = 0.706, P < 2e-16\n%d/%d concordant (99.7%%)", n_conc, n_de)) +
  labs(x = "log2 FC discovery (GSE106737)", y = "log2 FC replication\n(GSE83452)") +
  theme_ipf_lf(base = 7) +
  theme(plot.margin = margin(5, 8, 5, 5))

## ---------------- 组装：Row1 = a 设计 + b 火山（hero）；Row2 = c(宽)/d/e ----------------
banner_path <- "/mnt/results/fig1a_design_v2.png"
banner <- wrap_elements(full = cowplot::ggdraw() + cowplot::draw_image(banner_path))
design <- "
AAAAABBBBBBB
CCCCCDDDEEEE
"
combined <- banner + pB + pC + pD + pE +
  plot_layout(design = design, heights = c(0.85, 1.15)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "Liberation Sans"))

ggsave(file.path(out, "Fig1_main_v2.svg"), combined, width = 7.2, height = 4.4,
       device = svglite::svglite)
ggsave(file.path(out, "Fig1_main_v2.png"), combined, width = 7.2, height = 4.4, dpi = 600)

# 单面板导出（检查点用）
save_panel(pB, file.path(out, "fig1v2_b_volcano"), width = 3.6, height = 2.9)
save_panel(pC, file.path(out, "fig1v2_c_go"), width = 3.4, height = 2.9)
save_panel(pD, file.path(out, "fig1v2_d_slopes"), width = 3.4, height = 2.8)
save_panel(pE, file.path(out, "fig1v2_e_replication"), width = 3.6, height = 2.8)
cat("done\n")
