#!/usr/bin/env Rscript
# make_fig6.R — Figure 6 main (a-f) + supplements S1-S4: Regression Potential Score validation
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(cowplot)})
source("/workspace/fig1_theme.R")
OUT <- "/workspace/fig6/figs"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SC <- "/workspace/fig6/scores"; TB <- "/workspace/fig6/tables"
stat <- fread(file.path(TB, "fig6_validation_stats.csv"))
gv <- function(k) stat[check == k]
cls3 <- c(regressor = "#BFDA8E", stable = "#C8CACB", progressor = "#DCA7EB")
arm_cols <- c(vehicle = "#C8CACB", reversal = "#BFDA8E", sema = "#A7C0DF", lani = "#E2EA4D",
              chow = "#BEF3C5", prog = "#DCA7EB")

## ---------- panel a: schematic ----------
pa <- ggdraw() + draw_image("/mnt/results/fig6a_schematic_v2.png")

## ---------- panel b: GSE193066 paired histology ----------
pd <- fread(file.path(TB, "gse193066_paired_rps.csv"))
pd[, class := factor(class, levels = c("regressor", "stable", "progressor"))]
pd[, regressor := as.integer(class == "regressor")]
pb1 <- ggplot(pd, aes(class, dRPS, fill = class)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.75) +
  scale_fill_manual(values = cls3, guide = "none") +
  labs(x = NULL, y = "Delta RPS (visit 2 - 1)") + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5))
pb2 <- ggplot(pd, aes(dfib, dRPS, color = class)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_point(size = 1.8, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6, show.legend = FALSE) +
  scale_color_manual(values = cls3, name = NULL) +
  labs(x = "Delta fibrosis stage", y = NULL) + theme_ipf_lf() +
  theme(legend.position = "none")
pb <- wrap_elements(full = pb1 + pb2 + plot_layout(widths = c(1, 1.3)))

## ---------- panel c: mouse drug cohorts ----------
m1 <- fread(file.path(SC, "scores_GSE269493.csv"))
m1i <- m1[group %in% c("vehicle", "reversal", "sema", "lani", "chow")]
m1i[, group := factor(group, levels = c("chow", "vehicle", "reversal", "sema", "lani"))]
pc1 <- ggplot(m1i, aes(group, RPS, fill = group)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.1, alpha = 0.75) +
  scale_fill_manual(values = arm_cols, guide = "none") +
  labs(x = "GSE269493 (CDAA-HFD)", y = "RPS") + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
m2 <- fread(file.path(SC, "scores_GSE294630.csv"))
m2[, arm := fcase(treatment == "CDAA_Baseline_group_6_w", "baseline 6w",
                  treatment == "CDAA_vehicle_QD_6_12w", "vehicle 12w",
                  treatment == "CDAA_vehicle_QD_6_18w", "vehicle 18w",
                  treatment == "Sema_10_nmol_per_kg_QD_6_12w", "sema 12w",
                  treatment == "Sema_10_nmol_per_kg_QD_6_18w", "sema 18w",
                  treatment == "Chow_reversal_vehicle_QD_6_12w", "chow-rev 12w",
                  treatment == "Chow_reversal_vehicle_QD_6_18w", "chow-rev 18w",
                  treatment == "Chow_18w", "chow 18w")]
m2[, arm := factor(arm, levels = c("baseline 6w", "vehicle 12w", "vehicle 18w", "sema 12w",
                                   "sema 18w", "chow-rev 12w", "chow-rev 18w", "chow 18w"))]
c2cols <- c("baseline 6w" = "#DCA7EB", "vehicle 12w" = "#C8CACB", "vehicle 18w" = "#C8CACB",
            "sema 12w" = "#A7C0DF", "sema 18w" = "#A7C0DF",
            "chow-rev 12w" = "#BFDA8E", "chow-rev 18w" = "#BFDA8E", "chow 18w" = "#BEF3C5")
pc2 <- ggplot(m2, aes(arm, RPS, fill = arm)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.1, alpha = 0.75) +
  scale_fill_manual(values = c2cols, guide = "none") +
  labs(x = "GSE294630 (CDAA-HFD time course)", y = NULL) + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5))
pc <- wrap_elements(full = pc1 + pc2 + plot_layout(widths = c(1, 1.5)))

## ---------- panel d: GSE135251 staging (discovery) ----------
s <- fread(file.path(SC, "scores_GSE135251.csv"))
s[, stage := factor(paste0("F", fibrosis), levels = paste0("F", 0:4))]
pd1 <- ggplot(s, aes(stage, RPS)) +
  geom_boxplot(fill = "#A7C0C7", alpha = 0.75, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.4, color = "grey40") +
  labs(x = "Fibrosis stage (GSE135251)", y = "RPS") + theme_ipf_lf()
sm <- melt(s, id.vars = c("sample", "stage"), measure.vars = c("up", "down", "RPS"),
           variable.name = "arm", value.name = "score")
smu <- sm[, .(m = mean(score), se = sd(score) / sqrt(.N)), by = .(stage, arm)]
smu[, arm := factor(arm, levels = c("up", "down", "RPS"),
                    labels = c("Up arm (pro-regr.)", "Down arm (disease)", "RPS"))]
pd2 <- ggplot(smu, aes(as.numeric(stage), m, color = arm)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.12, linewidth = 0.4) +
  scale_color_manual(values = c("#BFDA8E", "#DCA7EB", "#A7C0C7"), name = NULL) +
  scale_x_continuous(breaks = 1:5, labels = paste0("F", 0:4)) +
  labs(x = NULL, y = "Mean score") + theme_ipf_lf() +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.5))
pd6 <- wrap_elements(full = pd1 + pd2 + plot_layout(widths = c(1, 1.3)))

## ---------- panel e: incremental value over baseline fibrosis ----------
roc_df <- rbindlist(lapply(1:nrow(pd), function(i) NULL))  # placeholder
fit0 <- glm(regressor ~ fib1, pd, family = binomial)
fit1 <- glm(regressor ~ fib1 + RPS1, pd, family = binomial)
library(pROC)
r0 <- roc(pd$regressor, predict(fit0, type = "response"), quiet = TRUE)
r1 <- roc(pd$regressor, predict(fit1, type = "response"), quiet = TRUE)
roc_plot <- data.frame(
  fpr = c(1 - r0$specificities, 1 - r1$specificities),
  tpr = c(r0$sensitivities, r1$sensitivities),
  model = rep(c(sprintf("Baseline fibrosis (AUC %.2f)", as.numeric(auc(r0))),
                sprintf("+ baseline RPS (AUC %.2f)", as.numeric(auc(r1)))),
              c(length(r0$sensitivities), length(r1$sensitivities))))
roc_plot$model <- factor(roc_plot$model, levels = c(grep("Baseline fibrosis", unique(roc_plot$model), value = TRUE),
                                                    grep("baseline RPS", unique(roc_plot$model), value = TRUE)))
pe <- ggplot(roc_plot, aes(fpr, tpr, color = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("#C8CACB", "#75A025"), name = NULL) +
  labs(x = "1 - specificity", y = "Sensitivity") + theme_ipf_lf() +
  theme(legend.position = c(0.98, 0.05), legend.justification = c(1, 0),
        legend.text = element_text(size = 7), legend.background = element_blank())

## ---------- panel f: external paired cohorts ----------
hp <- fread(file.path(TB, "gse70779_paired_rps.csv"))
hpl <- melt(hp, id.vars = "pair_id", measure.vars = c("PRE", "EOT"), variable.name = "tp", value.name = "RPS")
pf1 <- ggplot(hpl, aes(tp, RPS, group = pair_id)) +
  geom_line(alpha = 0.6, color = "grey50", linewidth = 0.4) +
  geom_point(aes(color = tp), size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c(PRE = "#DCA7EB", EOT = "#BFDA8E"), guide = "none") +
  labs(x = "GSE70779 (HCV cure)", y = "RPS") + theme_ipf_lf()
gp <- fread(file.path(TB, "gse83452_paired_rps.csv"))
gpl <- melt(gp, id.vars = "patient", measure.vars = c("RPS_base", "RPS_fu"), variable.name = "tp", value.name = "RPS")
gpl[, tp := factor(tp, levels = c("RPS_base", "RPS_fu"), labels = c("Baseline", "1 yr post-RYGB"))]
pf2 <- ggplot(gpl, aes(tp, RPS, group = patient)) +
  geom_line(alpha = 0.6, color = "grey50", linewidth = 0.4) +
  geom_point(aes(color = tp), size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c(Baseline = "#DCA7EB", `1 yr post-RYGB` = "#BFDA8E"), guide = "none") +
  labs(x = "GSE83452 (bariatric surgery)", y = NULL) + theme_ipf_lf() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
pf <- wrap_elements(full = pf1 + pf2 + plot_layout(widths = c(1, 1.2)))

fig6 <- (pa | pb) / (pc | pd6) / (pe | pf) + plot_annotation(tag_levels = "a")
save_panel(fig6, file.path(OUT, "fig6_main"), width = 13, height = 13.5)

## ---------- S1: score QC ----------
sg <- fread("/workspace/fig6/score_genes.csv")
s1a <- ggplot(sg[, .N, by = .(component, arm)], aes(component, N, fill = arm)) +
  geom_col(position = "stack", alpha = 0.75, width = 0.6) +
  scale_fill_manual(values = c(up = "#BFDA8E", down = "#DCA7EB"), name = "Arm") +
  labs(x = "Fig2 component", y = "Genes in strict core") + theme_ipf_lf()
cov <- data.table(cohort = c("GSE135251", "GSE70779", "GSE193066", "GSE269493", "GSE294630", "GSE83452"),
                  found = c(160, 160, 155, 160, 160, 160))
s1b <- ggplot(cov, aes(reorder(cohort, found), found, fill = found)) +
  geom_col(alpha = 0.75, width = 0.6) + coord_flip() +
  scale_fill_gradient(low = "#DCA7EB", high = "#BFDA8E", guide = "none") +
  labs(x = NULL, y = "Core genes detected (of 160 mapped)") + theme_ipf_lf()
s1 <- s1a + s1b + plot_annotation(tag_levels = "a")
save_panel(s1, file.path(OUT, "fig6S1_score_qc"), width = 9, height = 4)

## ---------- S2: GSE135251 component trajectories (Fig2 component sets) ----------
st <- fread("/mnt/results/Figure2/tables/staging_scores_GSE135251.csv")
stm <- melt(st, id.vars = c("gsm", "fibrosis"), measure.vars = c("reversal", "emergent"),
            variable.name = "component", value.name = "score")
stm[, stage := factor(paste0("F", fibrosis), levels = paste0("F", 0:4))]
s2 <- ggplot(stm, aes(stage, score, fill = component)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.6, position = position_dodge(0.7)) +
  scale_fill_manual(values = c(reversal = "#DCA7EB", emergent = "#BFDA8E"), name = "Fig2 component") +
  labs(x = "Fibrosis stage (GSE135251)", y = "ssGSEA score") + theme_ipf_lf() +
  theme(legend.position = "bottom")
save_panel(s2, file.path(OUT, "fig6S2_component_staging"), width = 7, height = 4.5)

## ---------- S3: mouse progression internal control ----------
m1p <- m1[group %in% c("chow", paste0("prog_", c(3, 6, 12, 20), "w"), "vehicle")]
m1p[, week := as.numeric(sub("prog_(\\d+)w", "\\1", group))]
m1p[group == "chow", week := 0]; m1p[group == "vehicle", week := 20]
m1p[, grp2 := fifelse(group == "chow", "chow", fifelse(group == "vehicle", "vehicle 20w", "CDAA-HFD"))]
s3 <- ggplot(m1p, aes(week, RPS, color = grp2)) +
  geom_jitter(width = 0.4, size = 1.6, alpha = 0.75) +
  geom_smooth(data = m1p[grp2 == "CDAA-HFD"], method = "lm", se = TRUE, color = "#DCA7EB", linewidth = 0.7) +
  scale_color_manual(values = c(chow = "#BEF3C5", `CDAA-HFD` = "#DCA7EB", `vehicle 20w` = "#C8CACB"), name = NULL) +
  labs(x = "Weeks on CDAA-HFD (GSE269493 progression)", y = "RPS") + theme_ipf_lf() +
  theme(legend.position = "bottom")
save_panel(s3, file.path(OUT, "fig6S3_progression_control"), width = 6.5, height = 4.5)

## ---------- S4: sensitivity (lenient core + z-variant) ----------
sen <- fread(file.path(SC, "scores_sensitivity.csv"))
senl <- melt(sen, id.vars = c("cohort", "sample", "strict"), measure.vars = c("lenient", "zvar"),
             variable.name = "variant", value.name = "score")
rho <- senl[, .(rho = cor(strict, score, method = "spearman")), by = .(cohort, variant)]
rho[, lab := sprintf("rho = %.3f", rho)]
s4 <- ggplot(senl, aes(strict, score)) +
  geom_point(size = 0.5, alpha = 0.4, color = "#A7C0C7") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_text(data = rho, aes(label = lab), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
            size = 2.8, family = "Liberation Sans", inherit.aes = FALSE) +
  facet_grid(variant ~ cohort, scales = "free",
             labeller = labeller(variant = c(lenient = "Lenient core (164)", zvar = "Signed mean-z variant"))) +
  labs(x = "RPS (strict core, ssGSEA)", y = "Alternative score") + theme_ipf_lf() +
  theme(strip.text = element_text(size = 7))
save_panel(s4, file.path(OUT, "fig6S4_sensitivity"), width = 9, height = 6)

cat("MAKE_FIG6_DONE\n")
