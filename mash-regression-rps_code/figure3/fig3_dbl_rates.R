#!/usr/bin/env Rscript
# Recompute per-sample doublet rates (human GSE298719 + mouse GSE305330) for S1
suppressMessages({ library(Seurat); library(scDblFinder); library(dplyr); library(Matrix) })
OUT <- "/workspace/fig3_out"

rate_for <- function(mat, px) {
  o <- CreateSeuratObject(mat, project = px, min.cells = 3, min.features = 200)
  pat <- if (grepl("FRPBiopsy", px)) "^MT-" else "^mt-"
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = pat)
  o <- subset(o, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 15)
  sce <- as.SingleCellExperiment(o)
  sce <- scDblFinder(sce, dbr.sd = 1, verbose = FALSE)
  data.frame(sample = px, cells_postqc = ncol(o),
             doublet_rate = mean(sce$scDblFinder.class == "doublet"))
}

read10x <- function(dir, px, feat_file = "_features.tsv.gz") {
  feat <- read.delim(gzfile(file.path(dir, paste0(px, feat_file))), header = FALSE)
  mat <- readMM(gzfile(file.path(dir, paste0(px, "_matrix.mtx.gz"))))
  rownames(mat) <- make.unique(feat[[2]])
  colnames(mat) <- paste0(px, "_c", seq_len(ncol(mat)))
  mat
}

res <- list()
# mouse
dM <- "/workspace/fig3_data/GSE305330"
for (px in sub("_features.tsv.gz$", "", list.files(dM, pattern = "_features.tsv.gz$")))
  res[[px]] <- rate_for(read10x(dM, px), px)
# human
dH <- "/workspace/fig3_data/GSE298719"
for (px in sub("_features.tsv.gz$", "", list.files(dH, pattern = "_features.tsv.gz$")))
  res[[px]] <- rate_for(read10x(dH, px), px)

rates <- bind_rows(res)
rates$species <- ifelse(grepl("FRPBiopsy", rates$sample), "Human", "Mouse")
print(rates)
write.csv(rates, file.path(OUT, "doublet_rates.csv"), row.names = FALSE)
cat("DONE\n")
