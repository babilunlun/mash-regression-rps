#!/usr/bin/env Rscript
# fig6_prep_score.R — Regression Potential Score (RPS): ssGSEA direction-concordance score
# RPS = ssGSEA(core genes UP in regression) - ssGSEA(core genes DOWN in regression)
# Core set = Fig2 strict conserved core (162); mouse via 1:1 orthologs, human direction signs.
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages({library(data.table); library(GSVA); library(edgeR)})
OUT <- "/workspace/fig6"; dir.create(file.path(OUT, "scores"), showWarnings = FALSE, recursive = TRUE)
T2 <- "/mnt/results/Figure2/tables"

## ---------- 1) score definition ----------
comp <- fread(file.path(T2, "components.csv"))
core <- readLines(file.path(T2, "core_genes_strict.txt"))
cdt <- comp[gene %in% core]
arm_up <- cdt[regr_logFC > 0]$gene      # up in regression = pro-regression arm
arm_down <- cdt[regr_logFC < 0]$gene    # down in regression = disease-program arm
cat("strict core:", length(core), "| arm_up:", length(arm_up), "| arm_down:", length(arm_down), "\n")
print(table(cdt$component))
fwrite(cdt[, .(gene, component, regr_logFC, arm = fifelse(regr_logFC > 0, "up", "down"))],
       file.path(OUT, "score_genes.csv"))

## mouse ortholog arms (human direction signs)
ort <- fread(file.path(T2, "orthologs_1to1.csv"))
m_up <- ort[human_symbol %in% arm_up]$symbol
m_down <- ort[human_symbol %in% arm_down]$symbol
cat("mouse ortholog arms: up", length(m_up), "down", length(m_down), "\n")

## ---------- 2) helpers ----------
dedupe <- function(mat, sym) {
  o <- order(rowMeans(mat), decreasing = TRUE)
  keep <- o[!duplicated(sym[o])]
  m <- mat[keep, , drop = FALSE]; rownames(m) <- sym[keep]; m
}
score_mat <- function(mat, up, down, prefix) {
  up2 <- intersect(up, rownames(mat)); down2 <- intersect(down, rownames(mat))
  cat(prefix, "genes found: up", length(up2), "/", length(up), "| down", length(down2), "/", length(down), "\n")
  par <- ssgseaParam(exprData = as.matrix(mat), geneSets = list(up = up2, down = down2))
  sc <- gsva(par, verbose = FALSE)
  data.table(sample = colnames(mat), up = sc["up", ], down = sc["down", ],
             RPS = sc["up", ] - sc["down", ])
}

## ---------- 3) human cohorts ----------
# --- GSE135251 (archived logCPM) ---
x <- readRDS("/mnt/results/Figure2/data/gse135251_processed.rds")
m135 <- dedupe(x$logcpm, x$symbol)
s135 <- score_mat(m135, arm_up, arm_down, "GSE135251")
s135 <- merge(s135, as.data.table(x$meta), by.x = "sample", by.y = "gsm")
fwrite(s135, file.path(OUT, "scores", "scores_GSE135251.csv"))

# --- GSE70779 HCV (archived log matrix, 9 pairs pre/EOT) ---
y <- readRDS("/mnt/results/Figure2/data/expr_gene_GSE70779.rds")
m707 <- dedupe(y, rownames(y))
s707 <- score_mat(m707, arm_up, arm_down, "GSE70779")
fwrite(s707, file.path(OUT, "scores", "scores_GSE70779.csv"))

# --- GSE193066 (GCT, author-normalized; log2(x+1)) ---
gct <- fread(file.path(OUT, "data", "GSE193066_HUn164.gct.gz"), skip = 2)
gn <- gct$gene_name
m193 <- as.matrix(gct[, -(1:2)]); m193 <- log2(m193 + 1)
m193 <- dedupe(m193, gn)
s193 <- score_mat(m193, arm_up, arm_down, "GSE193066")
meta193 <- fread(file.path(T2, "metadata_GSE193066.csv"))
s193 <- merge(s193, meta193, by.x = "sample", by.y = "title", all.x = TRUE)
fwrite(s193, file.path(OUT, "scores", "scores_GSE193066.csv"))

## ---------- 4) mouse cohorts (Ensembl IDs -> logCPM -> symbol) ---
load_counts <- function(f) {
  hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
  d <- fread(f, skip = 1, header = FALSE)
  if (length(hdr) == ncol(d) - 1) hdr <- c("id", hdr)  # GSE269493: first column unnamed
  setnames(d, hdr)
  ids <- d[[1]]; m <- as.matrix(d[, -1]); rownames(m) <- sub("\\..*$", "", ids); m
}
mouse_score <- function(counts, meta, prefix) {
  keep <- rowSums(counts) > 0
  lcpm <- cpm(counts[keep, ], log = TRUE, prior.count = 1)
  sym <- ort$symbol[match(rownames(lcpm), ort$ensembl)]
  ok <- !is.na(sym)
  m <- dedupe(lcpm[ok, ], sym[ok])
  s <- score_mat(m, m_up, m_down, prefix)
  s[, sample := as.character(sample)]
  meta2 <- copy(meta)[, col_id := as.character(col_id)]
  merge(s, meta2, by.x = "sample", by.y = "col_id")
}
meta493 <- fread(file.path(T2, "metadata_GSE269493.csv"))
c493 <- cbind(load_counts(file.path(OUT, "data", "GSE269493_intervention_counts.txt.gz")),
              load_counts(file.path(OUT, "data", "GSE269493_progression_counts.txt.gz")))
c493 <- c493[, !duplicated(colnames(c493))]
s493 <- mouse_score(c493, meta493, "GSE269493")
fwrite(s493, file.path(OUT, "scores", "scores_GSE269493.csv"))

meta630 <- fread(file.path(T2, "metadata_GSE294630.csv"))
c630 <- load_counts(file.path(OUT, "data", "GSE294630_processed.tsv.gz"))
s630 <- mouse_score(c630, meta630, "GSE294630")
fwrite(s630, file.path(OUT, "scores", "scores_GSE294630.csv"))

cat("PREP_SCORE_DONE\n")
