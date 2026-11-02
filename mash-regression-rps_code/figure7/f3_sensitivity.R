# f3_sensitivity.R — further F3-stratum validation from existing GSE193066 data
# (a) reproduce F3 AUC; (b) bootstrap 95% CI; (c) leave-one-out AUCs;
# (d) permutation P; (e) within-F3 dRPS~dfibrosis; (f) other baseline strata AUCs
suppressMessages({library(data.table)})
RES <- "/mnt/results/paired_cohort_validation"
set.seed(20261031)

d <- fread("/mnt/results/Figure6/tables/scores_GSE193066.csv")
pd <- d[order(visit), .(RPS1 = RPS[1], RPS2 = RPS[2], fib1 = fibrosis[1], fib2 = fibrosis[2],
                        nas1 = nas[1], nas2 = nas[2]), by = patient]
pd <- pd[!is.na(RPS2)]
pd[, `:=`(dRPS = RPS2 - RPS1, dfib = fib2 - fib1, dnas = nas2 - nas1)]
pd[, regr := as.integer(dfib < 0)]
cat("paired patients:", nrow(pd), "| regressors:", sum(pd$regr), "\n")

auc <- function(score, y) {  # rank-based AUC
  r <- rank(score, ties.method = "average")
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

out <- list()
add <- function(...) { r <- list(...); out[[length(out) + 1]] <<- r; invisible(NULL) }

## per-stratum AUCs (baseline RPS -> future regression)
for (f in sort(unique(pd$fib1))) {
  s <- pd[fib1 == f]
  if (length(unique(s$regr)) < 2 || sum(s$regr) < 2) {
    add(analysis = paste0("stratum F", f, " AUC"), n = nrow(s), regressors = sum(s$regr),
        estimate = NA, detail = "too few regressors")
    next
  }
  a <- auc(s$RPS1, s$regr)
  w <- suppressWarnings(wilcox.test(RPS1 ~ regr, s, exact = FALSE))
  add(analysis = paste0("stratum F", f, " AUC"), n = nrow(s), regressors = sum(s$regr),
      estimate = round(a, 4), detail = paste0("Wilcoxon P=", signif(w$p.value, 3)))
}

## F3 stratum focused analyses
f3 <- pd[fib1 == 3]
a3 <- auc(f3$RPS1, f3$regr)
cat("F3: n =", nrow(f3), " regressors =", sum(f3$regr), " AUC =", round(a3, 4), "\n")

## (b) bootstrap CI (stratified by outcome)
B <- 10000
boot <- replicate(B, {
  i1 <- sample(which(f3$regr == 1), replace = TRUE)
  i0 <- sample(which(f3$regr == 0), replace = TRUE)
  auc(c(f3$RPS1[i1], f3$RPS1[i0]), c(rep(1, length(i1)), rep(0, length(i0))))
})
ci <- quantile(boot, c(0.025, 0.975))
add(analysis = "F3 AUC bootstrap", n = nrow(f3), regressors = sum(f3$regr),
    estimate = round(a3, 4), detail = paste0("95% CI ", round(ci[1], 3), "-", round(ci[2], 3)))

## (c) leave-one-out AUCs
loo <- sapply(seq_len(nrow(f3)), function(i) auc(f3$RPS1[-i], f3$regr[-i]))
add(analysis = "F3 AUC leave-one-out", n = nrow(f3), regressors = sum(f3$regr),
    estimate = round(median(loo), 4),
    detail = paste0("LOO range ", round(min(loo), 3), "-", round(max(loo), 3)))

## (d) permutation P for AUC
perm <- replicate(B, auc(f3$RPS1, sample(f3$regr)))
p_perm <- (1 + sum(perm >= a3)) / (B + 1)
add(analysis = "F3 AUC permutation", n = nrow(f3), regressors = sum(f3$regr),
    estimate = round(a3, 4), detail = paste0("permutation P=", signif(p_perm, 3), " (B=1e4)"))

## (e) within-F3 longitudinal: dRPS ~ dfibrosis
sp <- suppressWarnings(cor.test(f3$dRPS, f3$dfib, method = "spearman", exact = FALSE))
add(analysis = "F3 dRPS~dfibrosis", n = nrow(f3), regressors = sum(f3$regr),
    estimate = round(unname(sp$estimate), 4), detail = paste0("Spearman P=", signif(sp$p.value, 3)))

## F3 regressors: dRPS > 0 ?
w <- suppressWarnings(wilcox.test(f3[regr == 1]$dRPS, mu = 0, exact = FALSE))
add(analysis = "F3 regressor dRPS>0", n = sum(f3$regr == 1), regressors = sum(f3$regr),
    estimate = round(median(f3[regr == 1]$dRPS), 4),
    detail = paste0("signed-rank P=", signif(w$p.value, 3), "; ",
                    sum(f3[regr == 1]$dRPS > 0), "/", sum(f3$regr == 1), " positive"))

st <- rbindlist(out, fill = TRUE)
fwrite(st, file.path(RES, "f3_sensitivity_stats.csv"))
print(st)
cat("F3_DONE\n")
