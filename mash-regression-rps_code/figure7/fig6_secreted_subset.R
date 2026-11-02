#!/usr/bin/env Rscript
# fig6_secreted_subset.R — Secreted-subset RPS: how much of the RPS signal is carried
# by the secreted component (the blood-translation feasibility question).
# Subset = 54 down-arm + 1 up-arm genes with secreted_class != "intracellular"
# (manuscript §2.8 secretome definition). Score = ssGSEA(up_sub) - ssGSEA(down_sub),
# identical construction to the full RPS. Null = 1000 random same-size down-arm subsets.
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages({library(data.table); library(GSVA); library(pROC)})
set.seed(2026)
TB <- "/mnt/results/Figure6/tables"
OUT <- "/workspace/fig6"; dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

## ---------- gene sets ----------
sg  <- fread(file.path(TB, "score_genes.csv"))                 # 160 RPS genes, arm labels
ann <- fread(file.path(TB, "rps160_secretome_annotation.csv")) # secretome classes
stopifnot(all(ann$gene %in% sg$gene))
sec  <- ann[secreted_class != "intracellular"]                 # 54 down + 1 up
strc <- ann[secreted_class == "secreted"]                      # 25 down + 1 up (strict)
cat("secretome set:", nrow(sec), "(down", sum(sec$arm == "down"), "/ up", sum(sec$arm == "up"), ")",
    "| strict:", nrow(strc), "\n")
up_sec   <- sec[arm == "up"]$gene
down_sec <- sec[arm == "down"]$gene
down_str <- strc[arm == "down"]$gene
down_all <- sg[arm == "down"]$gene                             # 138 down-arm (null pool)

## ---------- expression matrices (mirror fig6_prep_score.R) ----------
dedupe <- function(mat, sym) {
  o <- order(rowMeans(mat), decreasing = TRUE)
  keep <- o[!duplicated(sym[o])]
  m <- mat[keep, , drop = FALSE]; rownames(m) <- sym[keep]; m
}
# GSE193066 (GCT, author-normalized; log2(x+1))
gct <- fread(file.path(OUT, "data", "GSE193066_HUn164.gct.gz"), skip = 2)
m193 <- dedupe(log2(as.matrix(gct[, -(1:2)]) + 1), gct$gene_name)
# GSE135251 (archived logCPM)
x <- readRDS("/mnt/results/Figure2/data/gse135251_processed.rds")
m135 <- dedupe(x$logcpm, x$symbol)

## ---------- score one matrix with secreted + strict + 1000 null sets ----------
score_all <- function(mat, tag, B = 1000) {
  up_f   <- intersect(up_sec, rownames(mat))
  dsec_f <- intersect(down_sec, rownames(mat))
  dstr_f <- intersect(down_str, rownames(mat))
  pool   <- intersect(down_all, rownames(mat))
  cat(tag, "found: up", length(up_f), "| down_sec", length(dsec_f), "/", length(down_sec),
      "| down_strict", length(dstr_f), "/", length(down_str), "| null pool", length(pool), "\n")
  nulls <- replicate(B, sample(pool, length(dsec_f)), simplify = FALSE)
  gs <- c(list(up = up_f, sec = dsec_f, strict = dstr_f,
               full_up = intersect(sg[arm == "up"]$gene, rownames(mat)),
               full_down = pool), nulls)
  names(gs)[-(1:5)] <- paste0("null", seq_len(B))
  sc <- gsva(ssgseaParam(exprData = as.matrix(mat), geneSets = gs), verbose = FALSE)
  # subset RPS keeps the up-minus-down construction; nulls use the same fixed up set
  r <- sc["up", ]
  list(sec = r - sc["sec", ], strict = r - sc["strict", ],
       full = sc["full_up", ] - sc["full_down", ],
       null = sweep(sc[paste0("null", seq_len(B)), , drop = FALSE], 2, r,
                    function(x, y) y - x))  # up - null_i
}

s193 <- score_all(m193, "GSE193066")
s135 <- score_all(m135, "GSE135251")

## ---------- phenotypes ----------
ph193 <- fread(file.path(TB, "scores_GSE193066.csv"))[, .(sample, patient, visit, fibrosis, nas)]
ph135 <- fread(file.path(TB, "scores_GSE135251.csv"))[, .(sample, fibrosis, nas)]

res <- list()
add <- function(key, set, stat, p = NA, detail = "") res[[length(res) + 1]] <<-
  data.table(check = key, set = set, stat = stat, p = p, detail = detail)

pspearman <- function(x, y, z) {  # partial Spearman via residualized ranks (mirrors original)
  rx <- rank(x); ry <- rank(y); rz <- apply(as.matrix(z), 2, rank)
  ct <- cor.test(residuals(lm(rx ~ rz)), residuals(lm(ry ~ rz)), method = "pearson")
  c(rho = unname(ct$estimate), p = ct$p.value)
}

## ---------- GSE193066 paired histology ----------
build_pd <- function(scores) {
  d <- data.table(sample = names(scores), sub = as.numeric(scores))
  d <- merge(d, ph193, by = "sample")
  pd <- d[order(visit), .(sub1 = sub[1], sub2 = sub[2], fib1 = fibrosis[1], fib2 = fibrosis[2],
                          nas1 = nas[1], nas2 = nas[2]), by = patient]
  pd <- pd[!is.na(sub2)]
  pd[, `:=`(dsub = sub2 - sub1, dfib = fib2 - fib1, dnas = nas2 - nas1)]
  pd[, class := fcase(dfib < 0, "regressor", dfib > 0, "progressor", default = "stable")]
  pd
}
f3_stats <- function(pd) {
  f3 <- pd[fib1 == 3]; f3[, regr := as.integer(class == "regressor")]
  ci <- as.numeric(ci.auc(roc(f3$regr, f3$sub1, quiet = TRUE), method = "delong"))
  wp <- wilcox.test(sub1 ~ regr, data = f3)$p.value
  adjp <- coef(summary(glm(regr ~ sub1 + nas1, data = f3, family = binomial)))["sub1", "Pr(>|z|)"]
  list(auc = ci[2], lo = ci[1], hi = ci[3], wp = wp, adjp = adjp, n = nrow(f3))
}
paired_stats <- function(pd, tag, setname) {
  ct <- cor.test(pd$dsub, pd$dfib, method = "spearman", exact = FALSE)
  add(paste0(tag, "_dsub_vs_dfib"), setname, unname(ct$estimate), ct$p.value)
  pr <- pspearman(pd$dsub, pd$dfib, pd$dnas)
  add(paste0(tag, "_dsub_vs_dfib_partial_dnas"), setname, pr["rho"], pr["p"])
  ct2 <- cor.test(pd$dsub, pd$dnas, method = "spearman", exact = FALSE)
  add(paste0(tag, "_dsub_vs_dnas"), setname, unname(ct2$estimate), ct2$p.value)
  f <- f3_stats(pd)
  add(paste0(tag, "_F3_AUC"), setname, f$auc, f$wp,
      sprintf("95%%CI %.3f-%.3f; n=%d; NAS-adj P=%.4f", f$lo, f$hi, f$n, f$adjp))
  invisible(f$auc)
}
pd_sec <- build_pd(s193$sec); pd_str <- build_pd(s193$strict); pd_full <- build_pd(s193$full)
paired_stats(pd_sec, "b193", "secreted54")
paired_stats(pd_str, "b193", "strict26")
paired_stats(pd_full, "b193", "full160_anchor")   # must reproduce rho=-0.31 / AUC 0.754

## ---------- GSE135251 staging ----------
stage_rho <- function(scores) {
  d <- merge(data.table(sample = names(scores), sub = as.numeric(scores)), ph135, by = "sample")
  ct <- cor.test(d$sub, d$fibrosis, method = "spearman", exact = FALSE)
  c(rho = unname(ct$estimate), p = ct$p.value)
}
r135 <- stage_rho(s135$sec);  add("d135_sub_vs_stage", "secreted54", unname(r135["rho"]), unname(r135["p"]))
r135s <- stage_rho(s135$strict); add("d135_sub_vs_stage", "strict26", unname(r135s["rho"]), unname(r135s["p"]))
r135f <- stage_rho(s135$full);  add("d135_sub_vs_stage", "full160_anchor", unname(r135f["rho"]), unname(r135f["p"]))

# NAS is the real cross-sectional anchor (manuscript §2.8: full RPS~NAS rho=-0.30, P=8.6e-6)
nas_rho <- function(scores) {
  d <- merge(data.table(sample = names(scores), sub = as.numeric(scores)), ph135, by = "sample")
  ct <- cor.test(d$sub, d$nas, method = "spearman", exact = FALSE)
  pr <- pspearman(d$sub, d$nas, d$fibrosis)
  c(rho = unname(ct$estimate), p = as.numeric(ct$p.value),
    prho = unname(pr["rho"]), pp = as.numeric(pr["p"]))
}
n135 <- nas_rho(s135$sec)
add("d135_sub_vs_nas", "secreted54", n135["rho"], n135["p"],
    sprintf("partial|fibrosis rho=%.3f P=%.2e", n135["prho"], n135["pp"]))
n135s <- nas_rho(s135$strict)
add("d135_sub_vs_nas", "strict26", n135s["rho"], n135s["p"],
    sprintf("partial|fibrosis rho=%.3f P=%.2e", n135s["prho"], n135s["pp"]))
n135f <- nas_rho(s135$full)
add("d135_sub_vs_nas", "full160_anchor", n135f["rho"], n135f["p"],
    sprintf("partial|fibrosis rho=%.3f P=%.2e", n135f["prho"], n135f["pp"]))

## ---------- null distributions (random same-size down-arm subsets) ----------
null_dfib <- apply(s193$null, 1, function(v) {
  pd <- build_pd(v); suppressWarnings(cor(pd$dsub, pd$dfib, method = "spearman"))
})
null_f3 <- apply(s193$null, 1, function(v) {
  pd <- build_pd(v); f3 <- pd[fib1 == 3]
  as.numeric(auc(roc(as.integer(f3$class == "regressor"), f3$sub1, quiet = TRUE)))
})
null_135 <- apply(s135$null, 1, function(v) {
  d <- merge(data.table(sample = colnames(s135$null), sub = v), ph135, by = "sample")
  suppressWarnings(cor(d$sub, d$fibrosis, method = "spearman"))
})
null_135_nas <- apply(s135$null, 1, function(v) {
  d <- merge(data.table(sample = colnames(s135$null), sub = v), ph135, by = "sample")
  suppressWarnings(cor(d$sub, d$nas, method = "spearman"))
})
obs <- c(dfib = unname(cor(pd_sec$dsub, pd_sec$dfib, method = "spearman")),
         f3 = f3_stats(pd_sec)$auc, stage = unname(r135["rho"]), nas = unname(n135["rho"]))
pct <- c(dfib = mean(null_dfib <= obs["dfib"]),          # one-sided: as or more negative
         f3 = mean(null_f3 >= obs["f3"]),                # one-sided: as or larger AUC
         stage = mean(null_135 <= obs["stage"]),
         nas = mean(null_135_nas <= obs["nas"]))
add("null_percentile_dfib", "secreted54", pct["dfib"], NA,
    sprintf("obs rho=%.3f; null median %.3f [%.3f, %.3f]", obs["dfib"],
            median(null_dfib), quantile(null_dfib, 0.025), quantile(null_dfib, 0.975)))
add("null_percentile_f3auc", "secreted54", pct["f3"], NA,
    sprintf("obs AUC=%.3f; null median %.3f [%.3f, %.3f]", obs["f3"],
            median(null_f3), quantile(null_f3, 0.025), quantile(null_f3, 0.975)))
add("null_percentile_stage", "secreted54", pct["stage"], NA,
    sprintf("obs rho=%.3f; null median %.3f [%.3f, %.3f]", obs["stage"],
            median(null_135), quantile(null_135, 0.025), quantile(null_135, 0.975)))
add("null_percentile_nas", "secreted54", pct["nas"], NA,
    sprintf("obs rho=%.3f; null median %.3f [%.3f, %.3f]", obs["nas"],
            median(null_135_nas), quantile(null_135_nas, 0.025), quantile(null_135_nas, 0.975)))

stat <- rbindlist(res)
fwrite(stat, file.path(OUT, "tables", "rps_secreted_subset_stats.csv"))
fwrite(pd_sec, file.path(OUT, "tables", "rps_secreted_subset_paired.csv"))
fwrite(data.table(null_dfib = null_dfib, null_f3 = null_f3, null_stage = null_135,
                  null_nas = null_135_nas),
       file.path(OUT, "tables", "rps_secreted_subset_nulls.csv"))
print(stat[, .(check, set, stat, p)])
cat("\n--- details ---\n")
stat[detail != "", cat(sprintf("%s (%s): %s\n", check, set, detail), sep = "")]
cat("SUBSET_DONE\n")
