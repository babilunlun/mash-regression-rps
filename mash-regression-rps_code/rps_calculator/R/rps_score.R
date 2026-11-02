# rps_score.R — core scoring functions for the RPS online calculator
# Mirrors the manuscript pipeline (fig6_prep_score.R) exactly:
#   1) collapse duplicate gene symbols keeping the highest-expressed row
#   2) ssGSEA on the two core arms; RPS = ssGSEA(up arm) - ssGSEA(down arm)
#
# The ssGSEA implementation below is a pure-R replica of GSVA v2's algorithm
# (per-sample ranks, average ties truncated to integer, |rank|^0.25 weights,
# integrated ECDF difference, global score-range normalization), included so
# the app has no Bioconductor dependency. Validated against GSVA 2.0.0:
# max |difference| = 0 over all 328 arm scores of the 164-sample reference.
#
# Important: ssGSEA scores are divided by the global enrichment-score range
# computed across ALL samples in the input matrix, and genes constant across
# the input matrix are filtered. Scoring a sample alone is therefore
# degenerate. To keep uploaded samples exactly comparable to the published
# cohort scores, each upload is scored JOINTLY with the fixed GSE193066
# reference matrix (164 samples, 17,231-gene universe): the reference columns
# anchor both the gene filter and the normalizing range. Validated: all 164
# reference samples re-scored through this path match the published batch
# scores with max |difference| = 0.

## Load the 160-gene score definition (columns: gene, component, regr_logFC, arm)
rps_load_genes <- function(path = "data/score_genes.csv") {
  g <- read.csv(path, stringsAsFactors = FALSE)
  list(up = g$gene[g$arm == "up"],
       down = g$gene[g$arm == "down"],
       table = g)
}

## Collapse duplicate symbols, keeping the row with the highest mean expression
rps_dedupe <- function(mat, sym) {
  o <- order(rowMeans(mat), decreasing = TRUE)
  keep <- o[!duplicated(sym[o])]
  m <- mat[keep, , drop = FALSE]
  rownames(m) <- sym[keep]
  m
}

## Pure-R ssGSEA, exact replica of GSVA v2 (ssgseaParam/gsva) output
my_ssgsea <- function(mat, gene_sets, alpha = 0.25, normalize = TRUE) {
  mat <- as.matrix(mat)
  if (ncol(mat) > 1) {                       # filter constant genes (GSVA v2)
    v <- apply(mat, 1, function(x) length(unique(x)))
    mat <- mat[v > 1, , drop = FALSE]
  }
  n <- nrow(mat)
  R <- apply(mat, 2, function(x) as.integer(rank(x, ties.method = "average")))
  Ra <- abs(R)^alpha
  gs_idx <- lapply(gene_sets, function(S) {
    i <- match(S, rownames(mat)); i[!is.na(i)]
  })
  es <- sapply(seq_len(ncol(mat)), function(j) {
    geneRanking <- order(R[, j], decreasing = TRUE)
    sapply(gs_idx, function(idx) {
      pos <- match(idx, geneRanking)
      k <- length(pos)
      w <- Ra[geneRanking[pos], j]
      cdf_in <- sum(w * (n - pos + 1)) / sum(w)
      cdf_out <- (n * (n + 1) / 2 - sum(n - pos + 1)) / (n - k)
      cdf_in - cdf_out
    })
  })
  if (normalize) es <- es / (max(es) - min(es))
  rownames(es) <- names(gene_sets); colnames(es) <- colnames(mat)
  es
}

## Score a genes x samples matrix; returns data.frame(sample, up, down, RPS)
rps_score_matrix <- function(mat, up, down) {
  sc <- my_ssgsea(mat, list(up = intersect(up, rownames(mat)),
                            down = intersect(down, rownames(mat))))
  data.frame(sample = colnames(mat),
             up = sc["up", ], down = sc["down", ],
             RPS = sc["up", ] - sc["down", ],
             row.names = NULL)
}

## Parse a user-uploaded table into a named numeric vector (gene -> expression).
## Accepts any table with one gene-identifier column and one numeric column;
## column names are matched case-insensitively ("gene"/"symbol"/"expression").
rps_parse_input <- function(df) {
  nms <- tolower(names(df))
  gene_col <- which(nms %in% c("gene", "symbol", "gene_name", "genesymbol", "hgnc_symbol"))
  if (length(gene_col) == 0) gene_col <- which(vapply(df, is.character, logical(1)))
  if (length(gene_col) == 0) gene_col <- 1L
  gene_col <- gene_col[1]
  num_col <- which(vapply(df, is.numeric, logical(1)) & seq_along(df) != gene_col)
  if (length(num_col) == 0) stop("No numeric expression column found.")
  expr_col <- num_col[1]
  values <- as.numeric(df[[expr_col]])
  genes <- as.character(df[[gene_col]])
  ok <- !is.na(values) & !is.na(genes) & genes != ""
  list(values = values[ok], genes = genes[ok],
       gene_col = names(df)[gene_col], expr_col = names(df)[expr_col])
}

## Score one uploaded vector jointly with the fixed reference matrix.
## Duplicate symbols in the upload are collapsed keeping the highest value
## (equivalent to the highest-mean rule for a single column). Only genes
## present in the reference universe are used.
rps_score_upload <- function(values, genes, arms, ref_mat) {
  v <- tapply(values, genes, max)
  common <- intersect(rownames(ref_mat), names(v))
  up_mat <- matrix(v[common], ncol = 1,
                   dimnames = list(common, "uploaded_sample"))
  joint <- cbind(ref_mat[common, ], up_mat[common, , drop = FALSE])
  sc <- rps_score_matrix(joint, arms$up, arms$down)
  list(score = sc[sc$sample == "uploaded_sample", ],
       n_up = length(intersect(arms$up, common)),
       n_down = length(intersect(arms$down, common)),
       missing_up = setdiff(arms$up, common),
       missing_down = setdiff(arms$down, common),
       n_input_genes = length(v),
       n_universe = length(common))
}
