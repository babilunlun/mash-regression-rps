#!/usr/bin/env Rscript
# run_coloc.R — coloc.abf for MR-nominal gene x layer x outcome pairs (scoping FIXED)
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
library(data.table); library(coloc)

QC <- "/workspace/fig4/qc"; REG <- "/workspace/fig4/regions"; REXP <- "/workspace/fig4/regions_exp"

ometa <- list(
  finngen_fibrochir      = list(type = "cc", s = 2653/487866, N = 487866),
  finngen_nafld          = list(type = "cc", s = 3479/487866, N = 487866),
  gcst90091033_nafld_emerge = list(type = "cc", s = 8434/778614, N = 778614),
  gcst90267352_pdff      = list(type = "quant", N = 33235),
  gcst90590317_pdff      = list(type = "quant", N = 37589),
  chen_pred_masld        = list(type = "quant", N = 371405),
  chen_pred_pdff         = list(type = "quant", N = 36675)
)

load_exposure_region <- function(gene, layer) {
  g_ <- gene
  if (grepl("^INTERVAL", layer)) {
    acc <- sub("INTERVAL_pQTL:", "", layer)
    f <- file.path(REXP, paste0("INTERVAL_pQTL:", acc, "__", g_, ".tsv"))
    if (!file.exists(f)) return(NULL)
    d <- fread(f)
    d <- d[!is.na(beta) & !is.na(se) & se > 0]
    d[, chr := as.integer(chr)]
    list(d = d, type = "quant", N = 3301)
  } else if (layer == "eQTLGen_blood_eQTL") {
    f <- file.path(REXP, "eQTLGen_blood_eQTL__regions_hg38.tsv")
    if (!file.exists(f)) return(NULL)
    d0 <- fread(f)
    d <- d0[d0$gene == g_]
    if (nrow(d) < 20) return(NULL)
    d <- d[!is.na(panel_eaf) & panel_eaf > 0.001 & panel_eaf < 0.999]
    d[, `:=`(se = 1 / sqrt(2 * panel_eaf * (1 - panel_eaf) * (n + z^2)),
             beta = z / sqrt(2 * panel_eaf * (1 - panel_eaf) * (n + z^2)))]
    d <- d[!is.na(beta) & se > 0]
    list(d = d[, .(chr = chr38, pos = pos38, rsid, ea, oa, beta, se, p, eaf = panel_eaf)],
         type = "quant", N = 31684)
  } else NULL
}

run_coloc_pair <- function(gene, layer, outcome) {
  g_ <- gene
  ex <- load_exposure_region(g_, layer)
  if (is.null(ex)) return(NULL)
  oy0 <- fread(file.path(REG, paste0(outcome, ".tsv")))
  oy <- oy0[oy0$gene == g_]
  if (nrow(oy) < 20) return(NULL)
  pal <- c("A/T", "T/A", "C/G", "G/C")
  ex$d[, k := paste(chr, pos)]; oy[, k := paste(chr, pos)]
  m <- merge(ex$d, oy[, .(k, o_ea = ea, o_oa = oa, o_beta = beta, o_se = se, o_p = p, o_eaf = eaf)],
             by = "k")
  m[, pair := paste(ea, oa, sep = "/")]
  m <- m[!(pair %in% pal)]
  m[, flip := (ea == o_oa & oa == o_ea)]
  m <- m[(ea == o_ea & oa == o_oa) | flip]
  m[, by := fifelse(flip, -o_beta, o_beta)]
  m[, byse := o_se]
  m <- m[!is.na(beta) & !is.na(by) & se > 0 & byse > 0]
  if (nrow(m) < 30) return(NULL)
  om <- ometa[[outcome]]
  maf <- pmin(m$eaf, 1 - m$eaf)
  maf[is.na(maf)] <- 0.25
  d1 <- list(beta = m$beta, varbeta = m$se^2, N = ex$N, type = ex$type, snp = m$k, MAF = maf)
  if (ex$type == "quant") d1$sdY <- 1
  d2 <- list(beta = m$by, varbeta = m$byse^2, N = om$N, type = om$type, snp = m$k, MAF = maf)
  if (om$type == "cc") d2$s <- om$s else d2$sdY <- 1
  cc <- try(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5), silent = TRUE)
  if (inherits(cc, "try-error")) return(NULL)
  s <- cc$summary
  data.table(gene = gene, layer = layer, outcome = outcome, nsnps = s["nsnps"],
             PPH0 = s["PP.H0.abf"], PPH1 = s["PP.H1.abf"], PPH2 = s["PP.H2.abf"],
             PPH3 = s["PP.H3.abf"], PPH4 = s["PP.H4.abf"])
}

mr <- fread(file.path(QC, "mr_results.csv"))
todo <- unique(mr[p < 0.05, .(gene, layer, outcome)])
todo <- todo[grepl("INTERVAL|eQTLGen", layer)]
cat("coloc pairs to run:", nrow(todo), "\n")
res <- list()
for (i in seq_len(nrow(todo))) {
  r <- try(run_coloc_pair(todo$gene[i], todo$layer[i], todo$outcome[i]), silent = TRUE)
  if (!inherits(r, "try-error") && !is.null(r)) res[[length(res) + 1]] <- r
  if (i %% 25 == 0) cat(i, "/", nrow(todo), "\n")
}
coloc_res <- if (length(res)) rbindlist(res) else data.table()
fwrite(coloc_res, file.path(QC, "coloc_results.csv"))
print(coloc_res[order(-PPH4)])
cat("COLOC_DONE\n")
