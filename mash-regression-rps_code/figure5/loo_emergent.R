#!/usr/bin/env Rscript
# Leave-one-out emergent-program scores for the 3 target genes (all are inside
# the 172-gene mouse emergent set, so cor(target, emergent_prog) is circular).
.libPaths(c("/workspace/fig3_rlib/library", .libPaths()))
suppressPackageStartupMessages({library(Seurat); library(UCell)})

so <- readRDS("/workspace/fig3_out/mouse_sc_scored.rds")
df <- readRDS("/workspace/target3_mouse_df.rds")
genes <- c("Trp53inp1", "Slamf6", "Fn3k")
emg <- trimws(readLines("/workspace/fig3_out/genes_emergent_mouse.txt"))
emg <- emg[nzchar(emg)]
emg <- intersect(emg, rownames(so))
cat("emergent genes in object:", length(emg), "\n")

lyrs <- Layers(so, assay = "RNA")
cat("RNA layers (head):", paste(head(lyrs, 15), collapse = ", "), "\n")
if (length(grep("^data", lyrs)) > 1) {
  cat("joining split data layers...\n")
  so <- JoinLayers(so, assay = "RNA")
}
mat <- GetAssayData(so, assay = "RNA", layer = "data")
cat("matrix:", nrow(mat), "x", ncol(mat), "\n")

for (g in genes) {
  sig <- setdiff(emg, g)
  sc <- ScoreSignatures_UCell(mat, features = list(emergLO = sig))
  if (g == genes[1]) cat("UCell dims:", nrow(sc), "x", ncol(sc), "| cols:", paste(colnames(sc), collapse=","), "\n")
  df[[paste0("emergentLO_", g)]] <- as.numeric(sc[colnames(so), 1])
  cat(g, "LOO done; cor with full emergent_prog =",
      round(cor(df[[paste0("emergentLO_", g)]], df$emergent_prog, method = "spearman"), 4), "\n")
}
saveRDS(df, "/workspace/target3_mouse_df.rds")

myeloid_re <- "Kupffer|Monocyte|Dendritic|Neutrophil|SAMac|Macro"
samac_cells <- grepl("SAMac", df$celltype_fine)
strata <- list(all = rep(TRUE, nrow(df)),
               myeloid = grepl(myeloid_re, df$celltype) | samac_cells,
               SAMac = samac_cells,
               HSC = df$celltype == "HSC_myofibroblast")
cor_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 20 || sd(x) == 0 || sd(y) == 0) return(data.frame(rho = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(rho = unname(ct$estimate), p = ct$p.value, n = n)
}
grid3 <- expand.grid(stratum = names(strata), gene = genes, stringsAsFactors = FALSE)
res3 <- do.call(rbind, lapply(seq_len(nrow(grid3)), function(i) {
  g <- grid3$gene[i]; idx <- strata[[grid3$stratum[i]]]
  cbind(grid3[i, , drop = FALSE], score = "emergent_prog_LOO",
        cor_one(df[[g]][idx], df[[paste0("emergentLO_", g)]][idx]))
}))
res3$rho <- round(res3$rho, 3); res3$p <- signif(res3$p, 3)
cat("\n== LOO-corrected emergent correlations ==\n")
print(as.data.frame(res3), row.names = FALSE)
write.csv(res3, "/workspace/target3_mouse_celllevel_cor_emergentLO.csv", row.names = FALSE)
cat("DONE\n")
