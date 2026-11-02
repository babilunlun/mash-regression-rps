#!/usr/bin/env Rscript
# Human sc rebuild with checkpointing + memory hygiene (32 GB machine)
# Phase 1: read -> QC -> scDblFinder -> Normalize -> HVG -> save objs (slow, ~45 min)
# Phase 2: merge -> ScaleData -> PCA -> Harmony -> UMAP -> clusters -> save integrated
suppressMessages({
  library(Seurat); library(dplyr); library(Matrix); library(ggplot2)
})
DATA <- "/workspace/fig3_data/GSE298719"; OUT <- "/workspace/fig3_out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CK1 <- file.path(OUT, "human_objs_qc.rds")

if (!file.exists(CK1)) {
  feats <- list.files(DATA, pattern = "_features.tsv.gz$")
  prefixes <- sub("_features.tsv.gz$", "", feats)
  prefixes <- prefixes[grepl("FRPBiopsy", prefixes)]
  cat("human samples:", length(prefixes), "\n")

  read_sample <- function(px) {
    feat <- read.delim(gzfile(file.path(DATA, paste0(px, "_features.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)
    mat <- readMM(gzfile(file.path(DATA, paste0(px, "_matrix.mtx.gz"))))
    bc <- read.delim(gzfile(file.path(DATA, paste0(px, "_barcodes.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)[[1]]
    colnames(mat) <- bc
    rownames(mat) <- if (ncol(feat) >= 2) make.unique(feat[[2]]) else make.unique(feat[[1]])
    mat
  }
  mats <- lapply(prefixes, read_sample); names(mats) <- prefixes
  cat("cells per sample:\n"); print(sapply(mats, ncol))

  meta <- data.frame(sample = prefixes, stringsAsFactors = FALSE) %>%
    mutate(donor = sub("^(FRPBiopsy_[0-9]+)_.*$", "\\1", sample),
           sex = sub("^FRPBiopsy_[0-9]+_(Male|Female)_.*$", "\\1", sample),
           group = sub("^.*_(Healthy|MASLD)$", "\\1", sample))
  print(meta)

  objs <- lapply(seq_along(mats), function(i) {
    o <- CreateSeuratObject(mats[[i]], project = prefixes[i], min.cells = 3, min.features = 200)
    o$sample <- prefixes[i]; o$donor <- meta$donor[i]; o$sex <- meta$sex[i]; o$group <- meta$group[i]
    o
  })
  rm(mats); gc(verbose = FALSE)

  qc_stats <- lapply(objs, function(o) {
    o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
    data.frame(sample = o$sample[1], n_cells_raw = ncol(o),
               median_nFeature = median(o$nFeature_RNA), median_nCount = median(o$nCount_RNA),
               median_mt = median(o$percent.mt))
  }) |> bind_rows()
  print(qc_stats)
  write.csv(qc_stats, file.path(OUT, "qc_human_sc_raw.csv"), row.names = FALSE)

  objs <- lapply(objs, function(o) {
    o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
    subset(o, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 15)
  })
  suppressMessages(library(scDblFinder))
  objs <- lapply(objs, function(o) {
    sce <- as.SingleCellExperiment(o)
    sce <- scDblFinder(sce, dbr.sd = 1)
    o$doublet <- sce$scDblFinder.class
    o[, o$doublet == "singlet"]
  })
  cat("cells after QC+doublet:\n"); print(sapply(objs, ncol))

  objs <- lapply(objs, function(o) {
    o <- NormalizeData(o, verbose = FALSE)
    o <- FindVariableFeatures(o, nfeatures = 3000, verbose = FALSE)
    o
  })
  feats_int <- SelectIntegrationFeatures(objs, nfeatures = 3000)
  saveRDS(list(objs = objs, feats_int = feats_int), CK1)
  cat("PHASE1 DONE\n")
} else cat("PHASE1 checkpoint exists, skipping\n")

## ---- phase 2 ----
ck <- readRDS(CK1)
objs <- ck$objs; feats_int <- ck$feats_int; rm(ck); gc(verbose = FALSE)
m <- merge(objs[[1]], objs[-1])
rm(objs); gc(verbose = FALSE)
VariableFeatures(m) <- feats_int
m <- ScaleData(m, verbose = FALSE)
m <- RunPCA(m, npcs = 50, verbose = FALSE)
suppressMessages(library(harmony))
m <- RunHarmony(m, group.by.vars = "sample", verbose = FALSE)
m <- RunUMAP(m, reduction = "harmony", dims = 1:30, verbose = FALSE)
m <- FindNeighbors(m, reduction = "harmony", dims = 1:30, verbose = FALSE)
m <- FindClusters(m, resolution = 0.5, verbose = FALSE)
saveRDS(m, file.path(OUT, "human_sc_integrated.rds"))
cat("PHASE2 DONE; clusters:", length(levels(m$seurat_clusters)), "\n")
