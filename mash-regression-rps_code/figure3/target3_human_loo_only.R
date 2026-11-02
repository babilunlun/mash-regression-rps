#!/usr/bin/env Rscript
# LOO-only follow-up for human (a): the main script completed through correlations;
# this computes leave-one-out emergent scores for targets inside the emergent set.
.libPaths(c("/workspace/fig3_rlib/library", .libPaths()))
suppressPackageStartupMessages({library(Seurat); library(UCell)})

so <- readRDS("/mnt/shared-workspace/shared/human_rebuild/out/human_sc_scored.rds")
genes <- c("TP53INP1", "SLAMF6", "FN3K")
emg <- trimws(readLines("/workspace/fig3_out/genes_emergent_human.txt"))
emg <- emg[nzchar(emg)]
emg <- intersect(emg, rownames(so))
inside <- genes[genes %in% emg]
cat("emergent genes in object:", length(emg), "\n")
cat("targets inside emergent set:", ifelse(length(inside), paste(inside, collapse = ","), "none"), "\n")
if (!length(inside)) { cat("DONE (nothing to correct)\n"); quit(save = "no") }

m <- so@meta.data
lab_col <- grep("celltype_label", names(m), value = TRUE)[1]
DefaultAssay(so) <- "RNA"
lyrs <- Layers(so, assay = "RNA")
if (length(grep("^data", lyrs)) > 1) { cat("joining layers...\n"); so <- JoinLayers(so, assay = "RNA") }
mat <- GetAssayData(so, assay = "RNA", layer = "data")
expr <- as.data.frame(t(as.matrix(mat[inside, , drop = FALSE])))

myeloid_re <- "SAMac|Kupffer|Monocyte|Macro|cDC|pDC|DC|Neutro|Mast|myeloid"
samac <- grepl("SAMac", m[[lab_col]])
strata <- list(all = rep(TRUE, nrow(m)),
               myeloid = grepl(myeloid_re, m[[lab_col]]) | samac,
               SAMac = samac,
               HSC = grepl("HSC|Stellate|myofibro", m[[lab_col]], ignore.case = TRUE))
print(sapply(strata, sum))

cor_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 20 || sd(x) == 0 || sd(y) == 0) return(data.frame(rho = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(rho = unname(ct$estimate), p = ct$p.value, n = n)
}
res3 <- NULL
for (g in inside) {
  sc <- ScoreSignatures_UCell(mat, features = list(emergLO = setdiff(emg, g)))
  v <- as.numeric(sc[colnames(so), 1])
  for (sn in names(strata)) {
    idx <- strata[[sn]]
    res3 <- rbind(res3, data.frame(stratum = sn, gene = g, score = "emergent_prog_LOO",
                                   cor_one(expr[[g]][idx], v[idx])))
  }
  cat(g, "LOO done\n")
}
res3$rho <- round(res3$rho, 3); res3$p <- signif(res3$p, 3)
cat("\n== human LOO-corrected emergent correlations ==\n")
print(as.data.frame(res3), row.names = FALSE)
write.csv(res3, "/workspace/human_target3/target3_human_celllevel_cor_emergentLO.csv", row.names = FALSE)
cat("DONE\n")
