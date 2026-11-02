#!/usr/bin/env Rscript
## fig6_loo.R — RPS leave-one-out (TP53INP1 removed from down arm)
## Recompute scores for all 6 cohorts with FULL and LOO arms using identical pipeline,
## then recompute every validation statistic from fig6_analysis.R for both.
suppressMessages({library(data.table); library(GSVA); library(edgeR)})
OUT <- "/workspace/fig6_loo"; dir.create(OUT, showWarnings = FALSE)
dir.create(file.path(OUT, "scores"), showWarnings = FALSE)
T2 <- "/mnt/results/Figure2/tables"

## ---------- arms ----------
sg <- fread("/mnt/results/Figure6/tables/score_genes.csv")
arm_up   <- sg[arm == "up"]$gene
arm_down <- sg[arm == "down"]$gene
stopifnot("TP53INP1" %in% arm_down)
arm_down_loo <- setdiff(arm_down, "TP53INP1")
cat("arms full: up", length(arm_up), "down", length(arm_down),
    "| LOO down:", length(arm_down_loo), "\n")

ort <- fread(file.path(T2, "orthologs_1to1.csv"))
m_up   <- ort[human_symbol %in% arm_up]$symbol
m_down <- ort[human_symbol %in% arm_down]$symbol
m_down_loo <- ort[human_symbol %in% arm_down_loo]$symbol
cat("mouse arms full:", length(m_up), length(m_down), "| LOO down:", length(m_down_loo), "\n")

## ---------- helpers (identical to fig6_prep_score.R) ----------
dedupe <- function(mat, sym) {
  o <- order(rowMeans(mat), decreasing = TRUE)
  keep <- o[!duplicated(sym[o])]
  m <- mat[keep, , drop = FALSE]; rownames(m) <- sym[keep]; m
}
score_mat <- function(mat, up, down) {
  up2 <- intersect(up, rownames(mat)); down2 <- intersect(down, rownames(mat))
  par <- ssgseaParam(exprData = as.matrix(mat), geneSets = list(up = up2, down = down2))
  sc <- gsva(par, verbose = FALSE)
  list(scores = data.table(sample = colnames(mat), up = sc["up", ], down = sc["down", ],
                           RPS = sc["up", ] - sc["down", ]),
       n_up = length(up2), n_down = length(down2), mat = mat)
}
load_counts <- function(f) {
  hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
  d <- fread(f, skip = 1, header = FALSE)
  if (length(hdr) == ncol(d) - 1) hdr <- c("id", hdr)
  setnames(d, hdr)
  ids <- d[[1]]; m <- as.matrix(d[, -1]); rownames(m) <- sub("\\..*$", "", ids); m
}
mouse_mat <- function(counts) {
  keep <- rowSums(counts) > 0
  lcpm <- cpm(counts[keep, ], log = TRUE, prior.count = 1)
  sym <- ort$symbol[match(rownames(lcpm), ort$ensembl)]
  ok <- !is.na(sym)
  dedupe(lcpm[ok, ], sym[ok])
}

## ---------- load expression matrices (once) ----------
mats <- list()
x <- readRDS("/mnt/results/Figure2/data/gse135251_processed.rds")
mats$GSE135251 <- list(mat = dedupe(x$logcpm, x$symbol), sp = "human")
y <- readRDS("/mnt/results/Figure2/data/expr_gene_GSE70779.rds")
mats$GSE70779 <- list(mat = dedupe(y, rownames(y)), sp = "human")
gct <- fread("/workspace/fig6/data/GSE193066_HUn164.gct.gz", skip = 2)
gn <- gct$gene_name
m193 <- as.matrix(gct[, -(1:2)]); m193 <- log2(m193 + 1)
mats$GSE193066 <- list(mat = dedupe(m193, gn), sp = "human")
## GSE83452 via annotation db
suppressMessages(library(hugene20sttranscriptcluster.db))
f <- "/workspace/fig6/data/GSE83452_series_matrix.txt.gz"
ln <- readLines(f, n = 200); beg <- grep("!series_matrix_table_begin", ln)
d <- fread(f, skip = beg, header = TRUE)
map <- AnnotationDbi::select(hugene20sttranscriptcluster.db,
                             keys = as.character(d$ID_REF), keytype = "PROBEID", columns = "SYMBOL")
map <- as.data.table(map)[!is.na(SYMBOL)]; map <- map[!duplicated(PROBEID)]
d2 <- d[as.character(ID_REF) %in% map$PROBEID]
sym <- map$SYMBOL[match(as.character(d2$ID_REF), map$PROBEID)]
m <- as.matrix(d2[, -1]); if (max(m) > 100) m <- log2(m + 1)
mats$GSE83452 <- list(mat = dedupe(m, sym), sp = "human")
## mouse
c493 <- cbind(load_counts("/workspace/fig6/data/GSE269493_intervention_counts.txt.gz"),
              load_counts("/workspace/fig6/data/GSE269493_progression_counts.txt.gz"))
c493 <- c493[, !duplicated(colnames(c493))]
mats$GSE269493 <- list(mat = mouse_mat(c493), sp = "mouse")
c630 <- load_counts("/workspace/fig6/data/GSE294630_processed.tsv.gz")
mats$GSE294630 <- list(mat = mouse_mat(c630), sp = "mouse")

## ---------- score all cohorts x {full, LOO} ----------
S <- list()
for (nm in names(mats)) {
  mm <- mats[[nm]]
  if (mm$sp == "human") { up_f <- arm_up; dn_f <- arm_down; dn_l <- arm_down_loo
  } else { up_f <- m_up; dn_f <- m_down; dn_l <- m_down_loo }
  sf <- score_mat(mm$mat, up_f, dn_f)
  sl <- score_mat(mm$mat, up_f, dn_l)
  sf$scores[, sample := as.character(sample)]; sl$scores[, sample := as.character(sample)]
  S[[nm]] <- list(full = sf$scores, loo = sl$scores, mat = mm$mat,
                  n_down_full = sf$n_down, n_down_loo = sl$n_down)
  cat(nm, "scored. down-arm genes used:", sf$n_down, "->", sl$n_down, "\n")
  fwrite(sl$scores, file.path(OUT, "scores", paste0("scores_loo_", nm, ".csv")))
}

## sanity: recomputed full vs archived full
for (nm in names(S)) {
  arch <- fread(file.path("/mnt/results/Figure6/tables", paste0("scores_", nm, ".csv")))
  if (!"sample" %in% names(arch)) setnames(arch, 1, "sample")
  arch[, sample := as.character(sample)]
  mg <- merge(S[[nm]]$full, arch, by = "sample", suffixes = c("", "_arch"))
  cat(nm, "full-vs-archived RPS Spearman:", round(cor(mg$RPS, mg$RPS_arch, method = "spearman"), 5), "\n")
}

## ---------- metadata ----------
meta193 <- fread(file.path(T2, "metadata_GSE193066.csv"))
meta493 <- fread(file.path(T2, "metadata_GSE269493.csv"))[, col_id := as.character(col_id)]
meta630 <- fread(file.path(T2, "metadata_GSE294630.csv"))[, col_id := as.character(col_id)]
meta707 <- fread(file.path(T2, "metadata_GSE70779.csv"))
meta452 <- fread("/mnt/results/Figure1/tables/metadata_GSE83452_RYGB16.csv")
meta135 <- as.data.table(x$meta)

## ---------- validation stats (replicating fig6_analysis.R) ----------
valstats <- function(tag) {  # tag = "full" or "loo"
  res <- list()
  add <- function(key, ...) res[[key]] <<- data.table(check = key, ...)
  ## b) GSE193066
  d <- merge(S$GSE193066[[tag]], meta193, by.x = "sample", by.y = "title")
  pd <- d[order(visit), .(RPS1 = RPS[1], RPS2 = RPS[2], fib1 = fibrosis[1], fib2 = fibrosis[2],
                          nas1 = nas[1], nas2 = nas[2]), by = patient]
  pd <- pd[!is.na(RPS2)]
  pd[, `:=`(dRPS = RPS2 - RPS1, dfib = fib2 - fib1, dnas = nas2 - nas1)]
  pd[, class := fcase(dfib < 0, "regressor", dfib > 0, "progressor", default = "stable")]
  add("b_kw_class", stat = kruskal.test(dRPS ~ class, pd)$statistic,
      p = kruskal.test(dRPS ~ class, pd)$p.value,
      detail = paste(pd[, .N, by = class][order(class)]$N, collapse = "/"))
  add("b_regressor_vs_progressor",
      stat = wilcox.test(dRPS ~ class, pd[class %in% c("regressor", "progressor")])$statistic,
      p = wilcox.test(dRPS ~ class, pd[class %in% c("regressor", "progressor")])$p.value)
  add("b_regressor_vs_stable",
      stat = wilcox.test(dRPS ~ class, pd[class %in% c("regressor", "stable")])$statistic,
      p = wilcox.test(dRPS ~ class, pd[class %in% c("regressor", "stable")])$p.value)
  add("b_regressor_dRPS_gt0", stat = wilcox.test(pd[class == "regressor"]$dRPS, mu = 0)$statistic,
      p = wilcox.test(pd[class == "regressor"]$dRPS, mu = 0)$p.value)
  ct <- cor.test(pd$dRPS, pd$dfib, method = "spearman", exact = FALSE)
  add("b_dRPS_vs_dfib", stat = unname(ct$estimate), p = ct$p.value)
  ct2 <- cor.test(pd$dRPS, pd$dnas, method = "spearman", exact = FALSE)
  add("b_dRPS_vs_dnas", stat = unname(ct2$estimate), p = ct2$p.value)
  ## c) mouse
  m1 <- merge(S$GSE269493[[tag]], meta493, by.x = "sample", by.y = "col_id")
  for (arm in c("reversal", "sema", "lani")) {
    w <- wilcox.test(m1[group == arm]$RPS, m1[group == "vehicle"]$RPS)
    add(paste0("c_GSE269493_", arm, "_vs_vehicle"), stat = w$statistic, p = w$p.value,
        detail = paste0("median ", round(median(m1[group == arm]$RPS), 3), " vs ",
                        round(median(m1[group == "vehicle"]$RPS), 3)))
  }
  ct3 <- cor.test(m1[grepl("prog", group)]$RPS,
                  as.numeric(sub("prog_(\\d+)w", "\\1", m1[grepl("prog", group)]$group)),
                  method = "spearman", exact = FALSE)
  add("c_GSE269493_progression_decline", stat = unname(ct3$estimate), p = ct3$p.value)
  m2 <- merge(S$GSE294630[[tag]], meta630, by.x = "sample", by.y = "col_id")
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
  ## d) GSE135251 stage
  s <- merge(S$GSE135251[[tag]], meta135, by.x = "sample", by.y = "gsm")
  ct4 <- cor.test(s$RPS, s$fibrosis, method = "spearman", exact = FALSE)
  add("d_RPS_vs_stage", stat = unname(ct4$estimate), p = ct4$p.value)
  ct5 <- cor.test(s$up, s$fibrosis, method = "spearman", exact = FALSE)
  add("d_uparm_vs_stage", stat = unname(ct5$estimate), p = ct5$p.value)
  ct6 <- cor.test(s$down, s$fibrosis, method = "spearman", exact = FALSE)
  add("d_downarm_vs_stage", stat = unname(ct6$estimate), p = ct6$p.value)
  kw2 <- kruskal.test(RPS ~ factor(fibrosis), s)
  add("d_kw_stage", stat = kw2$statistic, p = kw2$p.value)
  ## e) incremental value
  pd[, regressor := as.integer(class == "regressor")]
  ct7 <- cor.test(pd$RPS1, pd$fib1, method = "spearman", exact = FALSE)
  add("e_orthogonality_RPS1_vs_fib1", stat = unname(ct7$estimate), p = ct7$p.value)
  fit0 <- glm(regressor ~ fib1, pd, family = binomial)
  fit1 <- glm(regressor ~ fib1 + RPS1, pd, family = binomial)
  lr <- anova(fit0, fit1, test = "LRT")
  add("e_LRT_add_RPS1", stat = lr$Deviance[2], p = lr$`Pr(>Chi)`[2])
  auc <- function(y, p) { r <- rank(p); (sum(r[y == 1]) - sum(y == 1) * (sum(y == 1) + 1) / 2) / (sum(y == 1) * sum(y == 0)) }
  set.seed(7); B <- 2000
  da <- replicate(B, {
    i <- sample(nrow(pd), replace = TRUE); b <- pd[i]
    if (length(unique(b$regressor)) < 2) return(NA_real_)
    f0 <- glm(regressor ~ fib1, b, family = binomial)
    f1 <- glm(regressor ~ fib1 + RPS1, b, family = binomial)
    auc(b$regressor, predict(f1, type = "response")) - auc(b$regressor, predict(f0, type = "response"))
  })
  add("e_AUC_fib1", stat = auc(pd$regressor, predict(fit0, type = "response")))
  add("e_AUC_fib1_plus_RPS1", stat = auc(pd$regressor, predict(fit1, type = "response")))
  add("e_deltaAUC_boot", stat = mean(da, na.rm = TRUE),
      detail = sprintf("95%%CI %.3f to %.3f", quantile(da, 0.025, na.rm = TRUE), quantile(da, 0.975, na.rm = TRUE)))
  ## f) external paired
  h <- merge(S$GSE70779[[tag]], meta707, by.x = "sample", by.y = "gsm")
  hp <- dcast(h, pair_id ~ timepoint, value.var = "RPS")
  hp[, dRPS := EOT - PRE]
  w2 <- wilcox.test(hp$EOT, hp$PRE, paired = TRUE)
  add("f_GSE70779_paired", stat = w2$statistic, p = w2$p.value,
      detail = paste0("median dRPS ", round(median(hp$dRPS), 3), "; ", sum(hp$dRPS > 0), "/", nrow(hp), " up"))
  g <- S$GSE83452[[tag]]
  gp <- merge(meta452[, .(patient, gsm_baseline, gsm_followup, response)], g[, .(sample, RPS)],
              by.x = "gsm_baseline", by.y = "sample")
  setnames(gp, "RPS", "RPS_base")
  gp <- merge(gp, g[, .(sample, RPS)], by.x = "gsm_followup", by.y = "sample")
  setnames(gp, "RPS", "RPS_fu")
  gp[, dRPS := RPS_fu - RPS_base]
  w3 <- wilcox.test(gp$RPS_fu, gp$RPS_base, paired = TRUE)
  add("f_GSE83452_paired", stat = w3$statistic, p = w3$p.value,
      detail = paste0("median dRPS ", round(median(gp$dRPS), 3), "; ", sum(gp$dRPS > 0), "/", nrow(gp), " up"))
  rbindlist(res, fill = TRUE)
}

st_full <- valstats("full"); st_loo <- valstats("loo")
cmp <- merge(st_full, st_loo, by = "check", suffixes = c("_full", "_loo"))
arch <- fread("/mnt/results/Figure6/tables/fig6_validation_stats.csv")
cmp <- merge(cmp, arch[, .(check, stat_arch = stat, p_arch = p)], by = "check", all.x = TRUE)
fwrite(cmp, file.path(OUT, "fig6_loo_validation_stats.csv"))
print(cmp[, .(check, stat_arch, stat_full, stat_loo, p_arch, p_full, p_loo)])

## ---------- target-level circularity: TP53INP1 expression vs RPS ----------
tg <- list()
for (nm in names(S)) {
  gene <- if (nm %in% c("GSE269493", "GSE294630")) "Trp53inp1" else "TP53INP1"
  mat <- S[[nm]]$mat
  if (!gene %in% rownames(mat)) { tg[[nm]] <- NULL; next }
  e <- mat[gene, colnames(mat)]
  sf <- S[[nm]]$full; sl <- S[[nm]]$loo
  tg[[nm]] <- data.table(cohort = nm,
    cor_full = cor(e, sf$RPS[match(colnames(mat), sf$sample)], method = "spearman"),
    cor_loo  = cor(e, sl$RPS[match(colnames(mat), sl$sample)], method = "spearman"))
}
tg <- rbindlist(tg)
fwrite(tg, file.path(OUT, "tp53inp1_vs_rps_loo.csv"))
print(tg)
cat("LOO_DONE\n")
