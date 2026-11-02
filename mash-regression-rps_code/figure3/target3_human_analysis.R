#!/usr/bin/env Rscript
# Human arm of (a): 3-target cell-type attribution on the rebuilt human sc object.
# Mirrors the mouse analysis. Run on human-rebuild (32 GB) after rebuild completes.
.libPaths(c("/workspace/fig3_rlib/library", .libPaths()))
suppressPackageStartupMessages({library(Seurat); library(dplyr); library(tidyr)})

OUT <- "/workspace/human_target3"
dir.create(OUT, showWarnings = FALSE)
genes <- c("TP53INP1", "SLAMF6", "FN3K")

so <- readRDS("/mnt/shared-workspace/shared/human_rebuild/out/human_sc_scored.rds")
m <- so@meta.data
cat("cells:", ncol(so), " genes:", nrow(so), "\n")
cat("metadata columns:\n"); print(names(m))

## ---- light anchor checks (rebuild script runs its own; double-check here) ----
lab_col <- grep("celltype_label", names(m), value=TRUE)[1]
if (is.na(lab_col)) lab_col <- grep("celltype", names(m), value=TRUE)[1]
cat("label column:", lab_col, "\n")
print(table(m[[lab_col]]))
grp_col <- grep("^group$|condition", names(m), value=TRUE)[1]
cat("group column:", grp_col, "\n"); print(table(m[[grp_col]]))
anchor_cells <- ncol(so) == 70424
samac_n <- sum(grepl("SAMac", m[[lab_col]]))
anchor_samac <- samac_n == 195
cat(sprintf("ANCHOR n_cells==70424: %s | SAMac==195 (got %d): %s\n",
            anchor_cells, samac_n, anchor_samac))
if (!anchor_cells || !anchor_samac) {
  cat("ANCHORS_FAILED - inspect label mapping before trusting output\n")
}

## ---- score columns ----
num_cols <- names(m)[sapply(m, is.numeric)]
rev_col <- grep("reversal", num_cols, value=TRUE)[1]
emg_col <- grep("emergent", num_cols, value=TRUE)[1]
sam_col <- grep("SAMac", num_cols, value=TRUE)[1]
hsc_col <- grep("HSC", num_cols, value=TRUE)[1]
cat("score columns:", rev_col, emg_col, sam_col, hsc_col, "\n")

## ---- expression ----
DefaultAssay(so) <- "RNA"
lyrs <- Layers(so, assay = "RNA")
if (length(grep("^data", lyrs)) > 1) { cat("joining layers...\n"); so <- JoinLayers(so, assay = "RNA") }
expr <- tryCatch(FetchData(so, vars = genes), error = function(e) NULL)
if (is.null(expr)) expr <- as.data.frame(t(as.matrix(GetAssayData(so, assay = "RNA", layer = "data")[genes, ])))
df <- bind_cols(as_tibble(m), as_tibble(expr))
cat("genes present:", paste(intersect(genes, names(df)), collapse = ", "), "\n")

## ---- HSC activation score if absent ----
if (is.na(hsc_col)) {
  suppressPackageStartupMessages(library(UCell))
  hsc_genes <- intersect(c("ACTA2","COL1A1","COL3A1","TIMP1","LOX","SPARC","PDGFRB","DES"), rownames(so))
  cat("computing HSC activation UCell from", length(hsc_genes), "markers\n")
  mat <- GetAssayData(so, assay = "RNA", layer = "data")
  sc <- ScoreSignatures_UCell(mat, features = list(HSCact = hsc_genes))
  df$HSCact <- as.numeric(sc[colnames(so), 1])
  hsc_col <- "HSCact"
}
scores <- c(rev_col, emg_col, sam_col, hsc_col)
scores <- scores[!is.na(scores)]
cat("scores used:", paste(scores, collapse = ", "), "\n")

## ---- strata ----
myeloid_re <- "SAMac|Kupffer|Monocyte|Macro|cDC|pDC|DC|Neutro|Mast|myeloid"
samac_cells <- grepl("SAMac", df[[lab_col]])
hsc_cells <- grepl("HSC|Stellate|myofibro", df[[lab_col]], ignore.case = TRUE)
strata <- list(all = rep(TRUE, nrow(df)),
               myeloid = grepl(myeloid_re, df[[lab_col]]) | samac_cells,
               SAMac = samac_cells,
               HSC = hsc_cells)
print(sapply(strata, sum))

## ---- per-celltype expression ----
ct_sum <- df %>% pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") %>%
  group_by(gene, .data[[lab_col]]) %>%
  summarise(n = n(), mean_expr = round(mean(expr), 4), pct_expr = round(100*mean(expr > 0), 1), .groups = "drop") %>%
  arrange(gene, desc(mean_expr))
names(ct_sum)[2] <- "celltype"
write.csv(ct_sum, file.path(OUT, "target3_human_expr_by_celltype.csv"), row.names = FALSE)

ct_grp <- df %>% pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") %>%
  group_by(gene, .data[[grp_col]], .data[[lab_col]]) %>%
  summarise(n = n(), mean_expr = round(mean(expr), 4), pct_expr = round(100*mean(expr > 0), 1), .groups = "drop")
names(ct_grp)[2:3] <- c("group", "celltype")
write.csv(ct_grp, file.path(OUT, "target3_human_expr_by_celltype_group.csv"), row.names = FALSE)

## ---- cell-level correlations ----
cor_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 20 || sd(x) == 0 || sd(y) == 0) return(data.frame(rho = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(rho = unname(ct$estimate), p = ct$p.value, n = n)
}
grid <- expand.grid(stratum = names(strata), gene = genes, score = scores, stringsAsFactors = FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  idx <- strata[[grid$stratum[i]]]
  cbind(grid[i, , drop = FALSE], cor_one(df[[grid$gene[i]]][idx], df[[grid$score[i]]][idx]))
}))
res$rho <- round(res$rho, 3); res$p <- signif(res$p, 3)
write.csv(res, file.path(OUT, "target3_human_celllevel_cor.csv"), row.names = FALSE)
print(as.data.frame(res), row.names = FALSE)

grid2 <- expand.grid(stratum = c("all", "myeloid"), group = unique(as.character(df[[grp_col]])),
                     gene = genes, score = scores, stringsAsFactors = FALSE)
res2 <- do.call(rbind, lapply(seq_len(nrow(grid2)), function(i) {
  idx <- strata[[grid2$stratum[i]]] & as.character(df[[grp_col]]) == grid2$group[i]
  cbind(grid2[i, , drop = FALSE], cor_one(df[[grid2$gene[i]]][idx], df[[grid2$score[i]]][idx]))
}))
res2$rho <- round(res2$rho, 3); res2$p <- signif(res2$p, 3)
write.csv(res2, file.path(OUT, "target3_human_celllevel_cor_by_group.csv"), row.names = FALSE)

## ---- LOO emergent (TP53INP1 is inside the 214-gene human emergent set) ----
emg_file <- "/workspace/fig3_out/genes_emergent_human.txt"
if (file.exists(emg_file)) {
  emg <- trimws(readLines(emg_file)); emg <- emg[nzchar(emg)]
  emg <- intersect(emg, rownames(so))
  cat("human emergent genes in object:", length(emg), "\n")
  inside <- genes[genes %in% emg]
  cat("targets inside emergent set:", ifelse(length(inside), paste(inside, collapse = ","), "none"), "\n")
  if (length(inside)) {
    suppressPackageStartupMessages(library(UCell))
    mat <- GetAssayData(so, assay = "RNA", layer = "data")
    for (g in inside) {
      sc <- ScoreSignatures_UCell(mat, features = list(emergLO = setdiff(emg, g)))
      df[[paste0("emergentLO_", g)]] <- as.numeric(sc[colnames(so), 1])
    }
    grid3 <- expand.grid(stratum = names(strata), gene = inside, stringsAsFactors = FALSE)
    res3 <- do.call(rbind, lapply(seq_len(nrow(grid3)), function(i) {
      g <- grid3$gene[i]; idx <- strata[[grid3$stratum[i]]]
      cbind(grid3[i, , drop = FALSE], score = "emergent_prog_LOO",
            cor_one(df[[g]][idx], df[[paste0("emergentLO_", g)]][idx]))
    }))
    res3$rho <- round(res3$rho, 3); res3$p <- signif(res3$p, 3)
    print(as.data.frame(res3), row.names = FALSE)
    write.csv(res3, file.path(OUT, "target3_human_celllevel_cor_emergentLO.csv"), row.names = FALSE)
  }
}
cat("DONE\n")
