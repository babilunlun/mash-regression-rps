#!/usr/bin/env Rscript
## figs_supp.R — supplementary figure candidates for the v2 added analyses
## Fig1-S4 stratified arms (a: RYGB vs lifestyle scatter; b: baseline RPS by response)
## Fig2-S8 persistent structure (prog_rho vs regr_logFC by class)
## Fig6-S5 RPS leave-one-out (a: full vs LOO -log10P; b: TP53INP1-RPS cor dumbbell)
suppressMessages({library(data.table); library(ggplot2)})
has_pw <- requireNamespace("patchwork", quietly = TRUE)
if (has_pw) library(patchwork) else library(gridExtra)
okabe <- c("#0072B2", "#009E73", "#E69F00", "#D55E00", "#999999", "#56B4E9")
base <- theme_bw(base_size = 9, base_family = "Liberation Sans") +
  theme(panel.grid.minor = element_blank(), legend.key.size = unit(3.5, "mm"),
        plot.margin = margin(4, 6, 4, 4))
OUT <- "/workspace/figs_supp"; dir.create(OUT, showWarnings = FALSE)
save_both <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".svg")), p, width = w, height = h, device = svglite::svglite)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 400)
}

## ---------------- Fig1-S4 ----------------
prog <- fread("/mnt/results/Figure1/tables/de_paired41.csv")
gcol <- "gene"
prog353 <- prog[prog$adj.P.Val < 0.05, ][[gcol]]
rygb <- fread("/mnt/results/Figure1/tables/de_stratified_rygb.csv")
life <- fread("/mnt/results/Figure1/tables/de_stratified_lifestyle.csv")
m <- merge(rygb[, .(gene, logFC_rygb = logFC)], life[, .(gene, logFC_life = logFC)], by = "gene")
m <- m[gene %in% prog353]
m[, dir := ifelse(gene %in% prog[adj.P.Val < 0.05 & logFC < 0][[gcol]], "Down in regression (250)", "Up in regression (103)")]
rho <- cor(m$logFC_rygb, m$logFC_life, method = "spearman")
sgn <- sum(sign(m$logFC_rygb) == sign(m$logFC_life))
p_a <- ggplot(m, aes(logFC_life, logFC_rygb, color = dir)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) + geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_point(size = 1.1, alpha = 0.75) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.6, size = 2.8, family = "Liberation Sans",
           label = sprintf("rho == %.2f*'; sign-concordant '*%d*'/'*353", rho, sgn), parse = TRUE) +
  scale_color_manual(values = c(okabe[1], okabe[4])) +
  labs(x = "Lifestyle arm logFC (post vs pre)", y = "RYGB arm logFC (post vs pre)", color = NULL) +
  base + theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
               legend.background = element_rect(fill = alpha("white", 0.7), color = NA))

br <- fread("/mnt/results/Figure1/tables/baseline_rps_response.csv")
br$reconstructed_response <- factor(br$reconstructed_response, levels = c("responder", "nonresponder"))
n_tab <- table(br$reconstructed_response)
wp <- wilcox.test(RPS ~ reconstructed_response, br)$p.value
p_b <- ggplot(br, aes(reconstructed_response, RPS, fill = reconstructed_response)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, linewidth = 0.35) +
  geom_jitter(width = 0.12, size = 1.1, alpha = 0.7) +
  annotate("text", x = 1.5, y = max(br$RPS) * 1.02, size = 2.8, family = "Liberation Sans",
           label = sprintf("Wilcoxon P = %.3f", wp)) +
  scale_fill_manual(values = c(okabe[2], okabe[3]), guide = "none") +
  scale_x_discrete(labels = c(sprintf("Responder\n(n=%d)", n_tab[1]), sprintf("Nonresponder\n(n=%d)", n_tab[2]))) +
  labs(x = NULL, y = "Baseline RPS") + base
save_both(p_a, "Fig1-S4-a_stratified_scatter", 3.4, 3.2)
save_both(p_b, "Fig1-S4-b_baseline_rps", 2.6, 3.2)
fig1s4 <- if (has_pw) { (p_a + p_b + plot_annotation(tag_levels = "a")) & theme(plot.tag = element_text(face = "bold", size = 10, family = "Liberation Sans")) } else { gridExtra::arrangeGrob(p_a, p_b, ncol = 2) }
save_both(fig1s4, "Fig1-S4_stratified_arms", 6.2, 3.2)

## ---------------- Fig2-S8 ----------------
comp <- fread("/mnt/results/Figure2/tables/components.csv")
pers <- comp[component == "persistent"]
pers[, regr_p := 2 * pt(-abs(regr_t), df = 40)]
pers[, class := fifelse(sign(prog_rho) != sign(regr_logFC),
                        fifelse(regr_p < 0.05, "resolving_trend", "arrested"),
                        fifelse(regr_p < 0.05, "still_drifting", "same_sign_ns"))]
cls_n <- pers[, .N, by = class]
lab <- setNames(trimws(sprintf("%s (n=%s)", c("Arrested", "Resolving trend", "Same-sign NS", "Still drifting"),
                        format(cls_n$N[match(c("arrested","resolving_trend","same_sign_ns","still_drifting"), cls_n$class)], big.mark = ",", trim = TRUE))),
                c("arrested","resolving_trend","same_sign_ns","still_drifting"))
pers$class <- factor(pers$class, levels = c("same_sign_ns","arrested","resolving_trend","still_drifting"))
p_c <- ggplot(pers[order(class)], aes(prog_rho, regr_logFC, color = class)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) + geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_point(size = 0.6, alpha = 0.45) +
  scale_color_manual(values = c("same_sign_ns" = "#999999", "arrested" = okabe[1],
                                "resolving_trend" = okabe[2], "still_drifting" = okabe[4]),
                     labels = lab, breaks = names(lab)) +
  guides(color = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  labs(x = "Progression association (Spearman rho, GSE135251)", y = "Regression logFC (GSE106737)", color = NULL) +
  base + theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 1),
               legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
save_both(p_c, "Fig2-S8_persistent_structure", 4.2, 3.4)
fwrite(pers[, .(gene, prog_rho, prog_fdr, regr_logFC, regr_t, regr_p, regr_fdr, class)],
       file.path(OUT, "persistent_structure_full.csv"))

## ---------------- Fig6-S5 ----------------
loo <- fread("/mnt/results/Figure6/tables/fig6_loo_validation_stats.csv")
lp <- loo[!is.na(p_full) & !is.na(p_loo)]
lp[, label := gsub("^[a-f]_", "", check)]
lp[, neglog_full := -log10(p_full)]; lp[, neglog_loo := -log10(p_loo)]
lim <- c(0, max(c(lp$neglog_full, lp$neglog_loo)) * 1.45)
p_d <- ggplot(lp, aes(neglog_full, neglog_loo)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_point(size = 1.4, color = okabe[1]) +
  geom_text(aes(label = label), size = 1.7, family = "Liberation Sans", hjust = -0.12, vjust = 0.5, check_overlap = TRUE) +
  coord_equal(xlim = lim, ylim = lim, clip = "off") +
  labs(x = "Full RPS: -log10(P)", y = "LOO RPS (TP53INP1 removed): -log10(P)") + base
save_both(p_d, "Fig6-S5-a_loo_pvalues", 3.6, 3.4)

tp <- fread("/mnt/results/Figure6/tables/tp53inp1_vs_rps_loo.csv")
tpl <- melt(tp, id.vars = "cohort", variable.name = "version", value.name = "cor")
tpl$version <- factor(tpl$version, levels = c("cor_full", "cor_loo"), labels = c("Full RPS", "LOO RPS"))
tpl$cohort <- factor(tpl$cohort, levels = tp$cohort[order(tp$cor_full)])
p_e <- ggplot(tpl, aes(cor, cohort)) +
  geom_line(aes(group = cohort), color = "grey70", linewidth = 0.5) +
  geom_point(aes(color = version), size = 1.8) +
  scale_color_manual(values = c("Full RPS" = "#999999", "LOO RPS" = okabe[4])) +
  labs(x = "TP53INP1 vs RPS (Spearman rho)", y = NULL, color = NULL) +
  base + theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 0),
               legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
save_both(p_e, "Fig6-S5-b_tp53inp1_dumbbell", 3.2, 3.0)
fig6s5 <- if (has_pw) { (p_d + p_e + plot_annotation(tag_levels = "a")) & theme(plot.tag = element_text(face = "bold", size = 10, family = "Liberation Sans")) } else { gridExtra::arrangeGrob(p_d, p_e, ncol = 2) }
save_both(fig6s5, "Fig6-S5_rps_loo", 7.0, 3.4)
cat("FIGS_DONE\n")
