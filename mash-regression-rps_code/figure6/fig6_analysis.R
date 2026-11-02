#!/usr/bin/env Rscript
# fig6_analysis.R — formal validation statistics for the Regression Potential Score (RPS)
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages(library(data.table))
OUT <- "/workspace/fig6"; SC <- file.path(OUT, "scores")
dir.create(file.path(OUT, "tables"), showWarnings = FALSE)
res <- list()
add <- function(key, ...) res[[key]] <<- data.table(check = key, ...)

## ---------- b) GSE193066 paired histology ----------
d <- fread(file.path(SC, "scores_GSE193066.csv"))
pd <- d[order(visit), .(RPS1 = RPS[1], RPS2 = RPS[2], up1 = up[1], up2 = up[2],
                        down1 = down[1], down2 = down[2],
                        fib1 = fibrosis[1], fib2 = fibrosis[2],
                        nas1 = nas[1], nas2 = nas[2]), by = patient]
pd <- pd[!is.na(RPS2)]
pd[, `:=`(dRPS = RPS2 - RPS1, dfib = fib2 - fib1, dnas = nas2 - nas1)]
pd[, class := fcase(dfib < 0, "regressor", dfib > 0, "progressor", default = "stable")]
fwrite(pd, file.path(OUT, "tables", "gse193066_paired_rps.csv"))

kw <- kruskal.test(dRPS ~ class, pd)
add("b_kw_class", stat = kw$statistic, p = kw$p.value,
    detail = paste(pd[, .N, by = class][order(class)]$N, collapse = "/"))
wr <- wilcox.test(dRPS ~ class, pd[class %in% c("regressor", "progressor")])
add("b_regressor_vs_progressor", stat = wr$statistic, p = wr$p.value)
ws <- wilcox.test(dRPS ~ class, pd[class %in% c("regressor", "stable")])
add("b_regressor_vs_stable", stat = ws$statistic, p = ws$p.value)
w1 <- wilcox.test(pd[class == "regressor"]$dRPS, mu = 0)
add("b_regressor_dRPS_gt0", stat = w1$statistic, p = w1$p.value)
ct <- cor.test(pd$dRPS, pd$dfib, method = "spearman", exact = FALSE)
add("b_dRPS_vs_dfib", stat = unname(ct$estimate), p = ct$p.value)
ct2 <- cor.test(pd$dRPS, pd$dnas, method = "spearman", exact = FALSE)
add("b_dRPS_vs_dnas", stat = unname(ct2$estimate), p = ct2$p.value)

## ---------- c) mouse drug cohorts ----------
m1 <- fread(file.path(SC, "scores_GSE269493.csv"))
for (arm in c("reversal", "sema", "lani")) {
  w <- wilcox.test(m1[group == arm]$RPS, m1[group == "vehicle"]$RPS)
  add(paste0("c_GSE269493_", arm, "_vs_vehicle"), stat = w$statistic, p = w$p.value,
      detail = paste0("median ", round(median(m1[group == arm]$RPS), 3), " vs ",
                      round(median(m1[group == "vehicle"]$RPS), 3)))
}
# progression internal control: RPS should decline 3w->20w
ct3 <- cor.test(m1[grepl("prog", group)]$RPS,
                as.numeric(sub("prog_(\\d+)w", "\\1", m1[grepl("prog", group)]$group)),
                method = "spearman", exact = FALSE)
add("c_GSE269493_progression_decline", stat = unname(ct3$estimate), p = ct3$p.value)

m2 <- fread(file.path(SC, "scores_GSE294630.csv"))
pairs630 <- list(c("Sema_10_nmol_per_kg_QD_6_12w", "CDAA_vehicle_QD_6_12w"),
                 c("Sema_10_nmol_per_kg_QD_6_18w", "CDAA_vehicle_QD_6_18w"),
                 c("Chow_reversal_vehicle_QD_6_12w", "CDAA_vehicle_QD_6_12w"),
                 c("Chow_reversal_vehicle_QD_6_18w", "CDAA_vehicle_QD_6_18w"))
for (pp in pairs630) {
  w <- wilcox.test(m2[treatment == pp[1]]$RPS, m2[treatment == pp[2]]$RPS)
  add(paste0("c_GSE294630_", sub("_QD.*", "", pp[1]), "_vs_", sub(".*QD_", "", pp[2])),
      stat = w$statistic, p = w$p.value,
      detail = paste0("median ", round(median(m2[treatment == pp[1]]$RPS), 3), " vs ",
                      round(median(m2[treatment == pp[2]]$RPS), 3)))
}

## ---------- d) GSE135251 staging (discovery) ----------
s <- fread(file.path(SC, "scores_GSE135251.csv"))
ct4 <- cor.test(s$RPS, s$fibrosis, method = "spearman", exact = FALSE)
add("d_RPS_vs_stage", stat = unname(ct4$estimate), p = ct4$p.value)
ct5 <- cor.test(s$up, s$fibrosis, method = "spearman", exact = FALSE)
add("d_uparm_vs_stage", stat = unname(ct5$estimate), p = ct5$p.value)
ct6 <- cor.test(s$down, s$fibrosis, method = "spearman", exact = FALSE)
add("d_downarm_vs_stage", stat = unname(ct6$estimate), p = ct6$p.value)
kw2 <- kruskal.test(RPS ~ factor(fibrosis), s)
add("d_kw_stage", stat = kw2$statistic, p = kw2$p.value)

## ---------- e) incremental value over baseline fibrosis (GSE193066) ----------
pd[, regressor := as.integer(class == "regressor")]
ct7 <- cor.test(pd$RPS1, pd$fib1, method = "spearman", exact = FALSE)
add("e_orthogonality_RPS1_vs_fib1", stat = unname(ct7$estimate), p = ct7$p.value)
fit0 <- glm(regressor ~ fib1, pd, family = binomial)
fit1 <- glm(regressor ~ fib1 + RPS1, pd, family = binomial)
lr <- anova(fit0, fit1, test = "LRT")
add("e_LRT_add_RPS1", stat = lr$Deviance[2], p = lr$`Pr(>Chi)`[2])
# AUCs with bootstrap CI for delta
auc <- function(y, p) {  # rank AUC
  r <- rank(p); (sum(r[y == 1]) - sum(y == 1) * (sum(y == 1) + 1) / 2) / (sum(y == 1) * sum(y == 0))
}
set.seed(7); B <- 2000
da <- replicate(B, {
  i <- sample(nrow(pd), replace = TRUE)
  b <- pd[i]
  if (length(unique(b$regressor)) < 2) return(NA_real_)
  f0 <- glm(regressor ~ fib1, b, family = binomial)
  f1 <- glm(regressor ~ fib1 + RPS1, b, family = binomial)
  auc(b$regressor, predict(f1, type = "response")) - auc(b$regressor, predict(f0, type = "response"))
})
add("e_AUC_fib1", stat = auc(pd$regressor, predict(fit0, type = "response")))
add("e_AUC_fib1_plus_RPS1", stat = auc(pd$regressor, predict(fit1, type = "response")))
add("e_deltaAUC_boot", stat = mean(da, na.rm = TRUE),
    detail = sprintf("95%%CI %.3f to %.3f", quantile(da, 0.025, na.rm = TRUE), quantile(da, 0.975, na.rm = TRUE)))

## ---------- f) external paired cohorts ----------
# GSE70779 HCV (9 pairs PRE -> EOT)
h <- fread(file.path(SC, "scores_GSE70779.csv"))
hm <- fread("/mnt/results/Figure2/tables/metadata_GSE70779.csv")
h <- merge(h, hm, by.x = "sample", by.y = "gsm")
hp <- dcast(h, pair_id ~ timepoint, value.var = "RPS")
hp[, dRPS := EOT - PRE]
w2 <- wilcox.test(hp$EOT, hp$PRE, paired = TRUE)
add("f_GSE70779_paired", stat = w2$statistic, p = w2$p.value,
    detail = paste0("median dRPS ", round(median(hp$dRPS), 3), "; ", sum(hp$dRPS > 0), "/", nrow(hp), " up"))
fwrite(hp, file.path(OUT, "tables", "gse70779_paired_rps.csv"))

# GSE83452 RYGB 16 pairs
g <- fread(file.path(SC, "scores_GSE83452.csv"))
gm <- fread("/mnt/results/Figure1/tables/metadata_GSE83452_RYGB16.csv")
gp <- merge(gm[, .(patient, gsm_baseline, gsm_followup, response)],
            g[, .(gsm, RPS)], by.x = "gsm_baseline", by.y = "gsm")
setnames(gp, "RPS", "RPS_base")
gp <- merge(gp, g[, .(gsm, RPS)], by.x = "gsm_followup", by.y = "gsm")
setnames(gp, "RPS", "RPS_fu")
gp[, dRPS := RPS_fu - RPS_base]
w3 <- wilcox.test(gp$RPS_fu, gp$RPS_base, paired = TRUE)
add("f_GSE83452_paired", stat = w3$statistic, p = w3$p.value,
    detail = paste0("median dRPS ", round(median(gp$dRPS), 3), "; ", sum(gp$dRPS > 0), "/", nrow(gp), " up"))
fwrite(gp, file.path(OUT, "tables", "gse83452_paired_rps.csv"))

stat <- rbindlist(res, fill = TRUE)
fwrite(stat, file.path(OUT, "tables", "fig6_validation_stats.csv"))
print(stat)
cat("ANALYSIS_DONE\n")
