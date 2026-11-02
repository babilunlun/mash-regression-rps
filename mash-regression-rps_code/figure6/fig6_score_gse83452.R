#!/usr/bin/env Rscript
# fig6_score_gse83452.R — score GSE83452 (RYGB 16 pairs) from series matrix
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages({library(data.table); library(GSVA)})
OUT <- "/workspace/fig6"

## arms (from prep script outputs)
sg <- fread(file.path(OUT, "score_genes.csv"))
arm_up <- sg[arm == "up"]$gene; arm_down <- sg[arm == "down"]$gene

## series matrix: data table starts after !series_matrix_table_begin
f <- file.path(OUT, "data", "GSE83452_series_matrix.txt.gz")
ln <- readLines(f, n = 200)
beg <- grep("!series_matrix_table_begin", ln)
d <- fread(f, skip = beg, header = TRUE)  # header line is ID_REF GSM...
cat("matrix:", nrow(d), "probes x", ncol(d) - 1, "samples\n")

## probe -> gene symbol
suppressMessages(library(hugene20sttranscriptcluster.db))
map <- AnnotationDbi::select(hugene20sttranscriptcluster.db,
                             keys = as.character(d$ID_REF), keytype = "PROBEID",
                             columns = "SYMBOL")
map <- as.data.table(map)[!is.na(SYMBOL)]
map <- map[!duplicated(PROBEID)]
d2 <- d[as.character(ID_REF) %in% map$PROBEID]
sym <- map$SYMBOL[match(as.character(d2$ID_REF), map$PROBEID)]
m <- as.matrix(d2[, -1])
cat("value range:", range(m[1:500, ]), "\n")  # check log scale
if (max(m) > 100) m <- log2(m + 1)
o <- order(rowMeans(m), decreasing = TRUE)
keep <- o[!duplicated(sym[o])]
m <- m[keep, ]; rownames(m) <- sym[keep]
cat("gene-level matrix:", nrow(m), "x", ncol(m), "\n")

## score
par <- ssgseaParam(exprData = m, geneSets = list(up = intersect(arm_up, rownames(m)),
                                                 down = intersect(arm_down, rownames(m))))
sc <- gsva(par, verbose = FALSE)
s <- data.table(gsm = colnames(m), up = sc["up", ], down = sc["down", ], RPS = sc["up", ] - sc["down", ])
cat("genes used: up", length(intersect(arm_up, rownames(m))), "down", length(intersect(arm_down, rownames(m))), "\n")
fwrite(s, file.path(OUT, "scores", "scores_GSE83452.csv"))
cat("GSE83452_DONE\n")
