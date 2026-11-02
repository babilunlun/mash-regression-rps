#!/usr/bin/env Rscript
# fig6_sensitivity.R — S4: RPS robustness — lenient core (164) + signed mean-z variant vs primary ssGSEA RPS
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages({library(data.table); library(GSVA); library(edgeR)})
OUT <- "/workspace/fig6"; T2 <- "/mnt/results/Figure2/tables"

comp <- fread(file.path(T2, "components.csv"))
corel <- readLines(file.path(T2, "core_genes_lenient.txt"))
ldt <- comp[gene %in% corel]
l_up <- ldt[regr_logFC > 0]$gene; l_down <- ldt[regr_logFC < 0]$gene
sg <- fread(file.path(OUT, "score_genes.csv"))
s_up <- sg[arm == "up"]$gene; s_down <- sg[arm == "down"]$gene
cat("lenient arms: up", length(l_up), "down", length(l_down), "\n")

dedupe <- function(mat, sym) {
  o <- order(rowMeans(mat), decreasing = TRUE)
  keep <- o[!duplicated(sym[o])]
  m <- mat[keep, , drop = FALSE]; rownames(m) <- sym[keep]; m
}
# signed mean-z variant: per-gene z across samples, oriented by regression direction
zscore_rps <- function(mat, dt) {
  g <- intersect(dt$gene, rownames(mat))
  z <- t(scale(t(mat[g, , drop = FALSE])))
  sgn <- dt$regr_logFC[match(g, dt$gene)] > 0
  colMeans(z[sgn, , drop = FALSE]) - colMeans(z[!sgn, , drop = FALSE])
}
ssgsea_rps <- function(mat, up, down) {
  par <- ssgseaParam(exprData = as.matrix(mat),
                     geneSets = list(up = intersect(up, rownames(mat)), down = intersect(down, rownames(mat))))
  sc <- gsva(par, verbose = FALSE)
  sc["up", ] - sc["down", ]
}

res <- list()
## human: GSE135251 + GSE193066
x <- readRDS("/mnt/results/Figure2/data/gse135251_processed.rds")
m135 <- dedupe(x$logcpm, x$symbol)
res[["GSE135251"]] <- data.table(sample = colnames(m135),
                                 strict = ssgsea_rps(m135, s_up, s_down),
                                 lenient = ssgsea_rps(m135, l_up, l_down),
                                 zvar = zscore_rps(m135, cdt <- comp[gene %in% readLines(file.path(T2, "core_genes_strict.txt"))]))

gct <- fread(file.path(OUT, "data", "GSE193066_HUn164.gct.gz"), skip = 2)
m193 <- dedupe(log2(as.matrix(gct[, -(1:2)]) + 1), gct$gene_name)
res[["GSE193066"]] <- data.table(sample = colnames(m193),
                                 strict = ssgsea_rps(m193, s_up, s_down),
                                 lenient = ssgsea_rps(m193, l_up, l_down),
                                 zvar = zscore_rps(m193, comp[gene %in% readLines(file.path(T2, "core_genes_strict.txt"))]))

## mouse: GSE269493 (ortholog-mapped)
ort <- fread(file.path(T2, "orthologs_1to1.csv"))
load_counts <- function(f) {
  hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
  d <- fread(f, skip = 1, header = FALSE)
  if (length(hdr) == ncol(d) - 1) hdr <- c("id", hdr)
  setnames(d, hdr)
  ids <- d[[1]]; m <- as.matrix(d[, -1]); rownames(m) <- sub("\\..*$", "", ids); m
}
c493 <- cbind(load_counts(file.path(OUT, "data", "GSE269493_intervention_counts.txt.gz")),
              load_counts(file.path(OUT, "data", "GSE269493_progression_counts.txt.gz")))
c493 <- c493[, !duplicated(colnames(c493))]
lcpm <- cpm(c493[rowSums(c493) > 0, ], log = TRUE, prior.count = 1)
sym <- ort$symbol[match(rownames(lcpm), ort$ensembl)]
ok <- !is.na(sym)
m493 <- dedupe(lcpm[ok, ], sym[ok])
ml_up <- ort[human_symbol %in% l_up]$symbol; ml_down <- ort[human_symbol %in% l_down]$symbol
ms_up <- ort[human_symbol %in% s_up]$symbol; ms_down <- ort[human_symbol %in% s_down]$symbol
# mouse z-variant: human core genes -> mouse symbols, human signs
cdt_m <- comp[gene %in% readLines(file.path(T2, "core_genes_strict.txt"))]
cdt_m <- merge(cdt_m, ort[, .(human_symbol, symbol)], by.x = "gene", by.y = "human_symbol")
cdt_m[, gene := symbol]
res[["GSE269493"]] <- data.table(sample = colnames(m493),
                                 strict = ssgsea_rps(m493, ms_up, ms_down),
                                 lenient = ssgsea_rps(m493, ml_up, ml_down),
                                 zvar = zscore_rps(m493, cdt_m))

rr <- rbindlist(res, idcol = "cohort")
fwrite(rr, file.path(OUT, "scores", "scores_sensitivity.csv"))
print(rr[, .(strict_vs_lenient = cor(strict, lenient, method = "spearman"),
             strict_vs_zvar = cor(strict, zvar, method = "spearman")), by = cohort])
cat("SENSITIVITY_DONE\n")
