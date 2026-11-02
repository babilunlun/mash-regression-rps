#!/usr/bin/env Rscript
# Fig2-S5: the MASH-specific layer of the regression signature
# = 353 MASH regression DE genes NOT direction-concordant in HCV cure (GSE70779)
# a) logFC scatter MASH vs HCV (concordant / non-concordant / ISG)
# b) component membership of non-concordant genes
# c) GO-BP enrichment of non-concordant genes (background = 353 MASH DE genes)
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages({ library(dplyr); library(ggplot2); library(patchwork); library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi) })
source("/workspace/fig1_theme.R")
OUT <- "/workspace/fig2s5"; dir.create(OUT, showWarnings = FALSE)

d <- read.csv("/mnt/results/Figure2/tables/hcv_vs_mash_regression.csv")
de <- d %>% filter(adj.P.Val_mash < 0.05)
cat("MASH DE genes:", nrow(de), "\n")
de <- de %>% mutate(concordant = sign(logFC_mash) == sign(logFC_hcv),
                    cls = case_when(is_isg ~ "ISG (antiviral)",
                                    concordant ~ "Shared (concordant in HCV)",
                                    TRUE ~ "MASH-specific (discordant in HCV)"))
print(table(de$cls))
cat("concordant excl ISG:", sum(de$concordant & !de$is_isg), "/", sum(!de$is_isg), "\n")

# component membership
comp <- read.csv("/mnt/results/Figure2/tables/components.csv") %>% dplyr::select(gene, component)
de <- de %>% left_join(comp, by = "gene")
de$component[is.na(de$component)] <- "not in universe"

# write the key table
de %>% filter(cls == "MASH-specific (discordant in HCV)") %>%
  arrange(adj.P.Val_mash) %>%
  write.csv(file.path(OUT, "mash_specific_regression_genes.csv"), row.names = FALSE)

## ---- a) scatter ----
set.seed(7)
cols <- c("Shared (concordant in HCV)" = "#C8CACB", "ISG (antiviral)" = "#B1C4D7",
          "MASH-specific (discordant in HCV)" = "#E2EA4D")
p_a <- ggplot(de, aes(x = logFC_mash, y = logFC_hcv)) +
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

## ---- b) component composition + persistent enrichment test ----
cmp <- de %>% count(component, cls) %>%
  group_by(cls) %>% mutate(frac = n / sum(n)) %>% ungroup()
cmp$component <- factor(cmp$component, levels = c("reversal", "emergent", "persistent", "other", "not in universe"))
# Fisher: persistent enrichment in MASH-specific vs shared
n_p_ms <- sum(de$cls == "MASH-specific (discordant in HCV)" & de$component == "persistent")
n_ms   <- sum(de$cls == "MASH-specific (discordant in HCV)")
n_p_sh <- sum(de$cls == "Shared (concordant in HCV)" & de$component == "persistent")
n_sh   <- sum(de$cls == "Shared (concordant in HCV)")
ft <- fisher.test(matrix(c(n_p_ms, n_ms - n_p_ms, n_p_sh, n_sh - n_p_sh), nrow = 2))
cat(sprintf("persistent: MASH-specific %d/%d (%.1f%%) vs shared %d/%d (%.1f%%), Fisher P=%.2e\n",
            n_p_ms, n_ms, 100 * n_p_ms / n_ms, n_p_sh, n_sh, 100 * n_p_sh / n_sh, ft$p.value))
cmp$cls <- factor(cmp$cls, levels = c("MASH-specific (discordant in HCV)", "Shared (concordant in HCV)", "ISG (antiviral)"))
p_b <- ggplot(cmp, aes(x = cls, y = n, fill = component)) +
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

## ---- c) GO enrichment: MASH-specific (non-ISG discordant) vs 353 background ----
bg <- de$gene
fg <- de %>% filter(cls == "MASH-specific (discordant in HCV)", !is_isg) %>% pull(gene)
cat("foreground (non-ISG MASH-specific):", length(fg), "\n")
bg_id <- bitr(bg, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
fg_id <- bitr(fg, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("mapped:", nrow(fg_id), "/", length(fg), "\n")
ego <- enrichGO(fg_id$ENTREZID, universe = bg_id$ENTREZID, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
eg <- as.data.frame(ego)
write.csv(eg, file.path(OUT, "mash_specific_go_bp.csv"), row.names = FALSE)
min_fdr <- min(eg$p.adjust)
cat("min GO FDR:", min_fdr, "\n")
# honest display: top nominal terms, labeled as non-FDR-significant if so
eg <- eg %>% arrange(pvalue) %>% head(10)
eg$Description <- factor(eg$Description, levels = rev(eg$Description))
fdr_note <- if (min_fdr >= 0.05) " (nominal; no term survives FDR)" else ""
p_c <- ggplot(eg, aes(x = -log10(pvalue), y = Description, fill = Count)) +
  geom_col(alpha = 0.85) +
  scale_fill_gradient(low = "#B1C4D7", high = "#5B7FA6") +
  theme_ipf_lf() + theme(legend.position = "bottom",
                         axis.text.y = element_text(size = 5.8),
                         plot.title = element_text(size = 7.5)) +
  labs(x = expression(-log[10](italic(P)[nominal])), y = NULL) +
  ggtitle(paste0("GO-BP, MASH-specific", fdr_note))

p <- p_a + p_b + p_c + plot_layout(widths = c(1.15, 0.75, 1.1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "Liberation Sans"))
ggsave(file.path(OUT, "Fig2-S5_mash_specific_layer.svg"), p, width = 11.5, height = 3.8, device = svglite::svglite)
ggsave(file.path(OUT, "Fig2-S5_mash_specific_layer.png"), p, width = 11.5, height = 3.8, dpi = 600)
cat("FIG2-S5 DONE\n")
