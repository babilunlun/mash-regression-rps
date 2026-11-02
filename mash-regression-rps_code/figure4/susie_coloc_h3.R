#!/usr/bin/env Rscript
# susie_coloc_h3.R — SuSiE-coloc re-examination of PPH3>0.8 ("LD-confounded") loci + ICAM1
# Multi-signal-aware: susie_rss per trait with EUR38s LD, then coloc.signals across signal pairs.
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table); library(coloc); library(susieR)
QC <- "/workspace/fig4/qc"; REG <- "/workspace/fig4/regions"; REXP <- "/workspace/fig4/regions_exp"
PANEL <- file.path(RAW <- "/workspace/fig4/raw", "EUR38s")
dir.create(file.path(QC, "ld"), showWarnings = FALSE)

ometa <- list(
  finngen_fibrochir = list(type = "cc", s = 2653/487866, N = 487866),
  finngen_nafld = list(type = "cc", s = 3479/487866, N = 487866),
  gcst90091033_nafld_emerge = list(type = "cc", s = 8434/778614, N = 778614),
  gcst90267352_pdff = list(type = "quant", N = 33235),
  gcst90590317_pdff = list(type = "quant", N = 37589),
  chen_pred_masld = list(type = "quant", N = 371405),
  chen_pred_pdff = list(type = "quant", N = 36675))

mr <- fread(file.path(QC, "mr_results.csv"))
cc <- fread(file.path(QC, "coloc_results.csv"))
todo <- merge(mr[fdr < 0.05, .(gene, layer, outcome)], cc[PPH3 > 0.8, .(gene, layer, outcome)],
              by = c("gene", "layer", "outcome"))
todo <- rbind(todo, data.table(gene = "ICAM1", layer = "eQTLGen_blood_eQTL", outcome = "chen_pred_masld"))
todo <- unique(todo)
cat("loci to re-examine:", nrow(todo), "\n")

get_ld <- function(snps, tag) {
  tf <- tempfile(tmpdir = file.path(QC, "ld"))
  fwrite(data.table(snps), paste0(tf, ".snps"), col.names = FALSE)
  system2("plink2", c("--pfile", PANEL, "--extract", paste0(tf, ".snps"),
                      "--export", "A-transpose", "--out", tf, "--silent"), stdout = FALSE, stderr = FALSE)
  tr <- fread(paste0(tf, ".traw"))
  unlink(paste0(tf, "*"))
  g <- as.matrix(tr[, -(1:6), with = FALSE])
  rownames(g) <- tr$SNP
  g <- g[snps, , drop = FALSE]
  # mean-impute missing, standardize
  for (i in seq_len(nrow(g))) {
    mu <- mean(g[i, ], na.rm = TRUE); g[i, is.na(g[i, ])] <- mu
  }
  sdv <- apply(g, 1, sd)
  keep <- sdv > 0
  LD <- cor(t(g[keep, , drop = FALSE]))
  # LDshrink-style shrinkage for reference-panel LD + nearest PD projection
  w <- 2 / ncol(g)
  LDs <- (1 - w) * LD + w * diag(nrow(LD))
  LDs <- as.matrix(Matrix::nearPD(LDs, corr = TRUE)$mat)
  dimnames(LDs) <- dimnames(LD)
  list(ld = LDs, kept = rownames(g)[keep])
}

# global lookups (read once)
PVAR <- fread(file.path(RAW, "EUR38s.pvar"), skip = "#CHROM")
setnames(PVAR, c("chr", "pos", "id", "ref", "alt"))
PVAR <- PVAR[chr %in% c(as.character(1:22), "X")]
PVAR[, chr := as.integer(fifelse(chr == "X", "23", chr))]
setkey(PVAR, chr, pos)
EQREG <- fread(file.path(REXP, "eQTLGen_blood_eQTL__regions_hg38.tsv"))

run_susie_coloc <- function(gene, layer, outcome) {
  g_ <- gene
  # exposure region
  if (grepl("^INTERVAL", layer)) {
    acc <- sub("INTERVAL_pQTL:", "", layer)
    f <- file.path(REXP, paste0("INTERVAL_pQTL:", acc, "__", g_, ".tsv"))
    if (!file.exists(f)) return(NULL)
    ex <- fread(f); ex[, chr := as.integer(chr)]
    exN <- 3301; extype <- "quant"
  } else {
    ex <- EQREG[EQREG$gene == g_]
    ex <- ex[!is.na(panel_eaf) & panel_eaf > 0.001 & panel_eaf < 0.999]
    ex[, `:=`(se = 1 / sqrt(2 * panel_eaf * (1 - panel_eaf) * (n + z^2)),
              beta = z / sqrt(2 * panel_eaf * (1 - panel_eaf) * (n + z^2)))]
    ex <- ex[, .(chr = chr38, pos = pos38, ea, oa, beta, se, eaf = panel_eaf)]
    exN <- 31684; extype <- "quant"
  }
  oy0 <- fread(file.path(REG, paste0(outcome, ".tsv")))
  oy <- oy0[oy0$gene == g_]
  ex[, k := paste(chr, pos)]; oy[, k := paste(chr, pos)]
  m <- merge(ex, oy[, .(k, o_ea = ea, o_oa = oa, o_beta = beta, o_se = se)], by = "k")
  pal <- c("A/T", "T/A", "C/G", "G/C")
  m <- m[!(paste(ea, oa, sep = "/") %in% pal)]
  m[, flip := (ea == o_oa & oa == o_ea)]
  m <- m[(ea == o_ea & oa == o_oa) | flip]
  m[, by := fifelse(flip, -o_beta, o_beta)]
  m <- m[!is.na(beta) & !is.na(by) & se > 0 & o_se > 0]
  if (nrow(m) < 50) return(NULL)
  # LD from panel: match by chr:pos
  m[, `:=`(c = as.integer(sub(" .*", "", k)), p = as.integer(sub(".* ", "", k)))]
  m[, pid := PVAR[.(m$c, m$p), id]]
  m <- m[!is.na(pid)]
  m <- m[!duplicated(pid)]
  if (nrow(m) < 50) return(NULL)
  ldres <- get_ld(m$pid, paste0(g_, "_", outcome))
  m <- m[pid %in% ldres$kept]
  ord <- match(ldres$kept, m$pid)
  m <- m[ord]
  LD <- ldres$ld
  om <- ometa[[outcome]]
  d1 <- list(beta = m$beta, varbeta = m$se^2, N = exN, type = extype, snp = m$pid, LD = LD)
  if (extype == "quant") d1$sdY <- 1
  d2 <- list(beta = m$by, varbeta = m$o_se^2, N = om$N, type = om$type, snp = m$pid, LD = LD)
  if (om$type == "cc") d2$s <- om$s else d2$sdY <- 1
  s1 <- try(runsusie(d1), silent = TRUE)
  s2 <- try(runsusie(d2), silent = TRUE)
  erv <- FALSE
  if (inherits(s1, "try-error") || inherits(s2, "try-error")) {
    # LD/summary-stat inconsistency fallback: estimate residual variance
    erv <- TRUE
    s1 <- try(runsusie(d1, estimate_residual_variance = TRUE), silent = TRUE)
    s2 <- try(runsusie(d2, estimate_residual_variance = TRUE), silent = TRUE)
  }
  if (inherits(s1, "try-error") || inherits(s2, "try-error")) return(NULL)
  if (is.null(s1$sets$cs) || is.null(s2$sets$cs)) {
    return(data.table(gene = gene, layer = layer, outcome = outcome,
                      n_cs_exposure = length(s1$sets$cs), n_cs_outcome = length(s2$sets$cs),
                      max_PPH4 = NA_real_, erv = erv, verdict = "no_CS_in_one_trait"))
  }
  res <- try(coloc.susie(s1, s2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5), silent = TRUE)
  if (inherits(res, "try-error")) return(NULL)
  sm <- res$summary
  best <- sm[which.max(sm$PP.H4.abf)]
  data.table(gene = gene, layer = layer, outcome = outcome,
             n_cs_exposure = length(s1$sets$cs), n_cs_outcome = length(s2$sets$cs),
             max_PPH4 = best$PP.H4.abf, erv = erv,
             verdict = ifelse(best$PP.H4.abf > 0.5, "colocalised_signal_pair", "distinct_signals"))
}

res <- list()
for (i in seq_len(nrow(todo))) {
  r <- try(run_susie_coloc(todo$gene[i], todo$layer[i], todo$outcome[i]), silent = TRUE)
  if (!inherits(r, "try-error") && !is.null(r)) res[[length(res) + 1]] <- r
  cat(i, "/", nrow(todo), todo$gene[i], todo$outcome[i],
      ifelse(!inherits(r, "try-error") && !is.null(r), paste("PPH4 =", signif(r$max_PPH4, 2), r$verdict), "failed"), "\n")
}
out <- rbindlist(res, fill = TRUE)
fwrite(out, file.path(QC, "susie_coloc_h3_recheck.csv"))
print(out)
cat("SUSIE_COLOC_DONE\n")
