#!/usr/bin/env Rscript
# Figure 3 - mouse scRNA (GSE305330): QC, integration, clustering, cell-type annotation
# Design: 2 regression models (CCl4 withdrawal, FAT-MASH dietary) x 3 timepoints
# (peak, 2wk regression, 6wk regression). CellBender-filtered 10X matrices.
suppressMessages({
  library(Seurat); library(dplyr); library(Matrix); library(ggplot2)
})
DATA <- "/workspace/fig3_data/GSE305330"; OUT <- "/workspace/fig3_out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- read samples ----
# per-sample files: <prefix>_features.tsv.gz, <prefix>_matrix.mtx.gz,
# <prefix>_cellbender_filtered_feature_bc_matrix_cell_barcodes.csv.gz (and possibly _barcodes.tsv.gz)
feats <- list.files(DATA, pattern = "_features.tsv.gz$")
prefixes <- sub("_features.tsv.gz$", "", feats)
cat("samples found:", length(prefixes), "\n")

read_sample <- function(px) {
  f <- file.path(DATA, paste0(px, "_features.tsv.gz"))
  m <- file.path(DATA, paste0(px, "_matrix.mtx.gz"))
  feat <- read.delim(gzfile(f), header = FALSE, stringsAsFactors = FALSE)
  mat <- readMM(gzfile(m))
  # matrix is already CellBender-filtered; no matching barcode file uploaded -> synthetic names
  colnames(mat) <- paste0(px, "_cell", seq_len(ncol(mat)))
  rownames(mat) <- if (ncol(feat) >= 2) make.unique(feat[[2]]) else make.unique(feat[[1]])
  mat
}

mats <- lapply(prefixes, read_sample)
names(mats) <- prefixes
cat("cells per sample:\n"); print(sapply(mats, ncol))

# condition metadata from sample names: KL70-CCl6wk / KL70-FAT8wk etc are PROGRESSION naming;
# series matrix titles give the true labels (Peak/2wk/6wk Regression). Map by GSM order.
sm <- readLines(gzfile("/workspace/fig3_data/GSE305330_series_matrix.txt.gz"))
srow <- sm[grepl("^!Sample_title", sm)]
titles <- gsub('"', "", strsplit(srow, "\t")[[1]][-1])
grow <- sm[grepl("^!Sample_geo_accession", sm)]
gsms <- gsub('"', "", strsplit(grow, "\t")[[1]][-1])
cat("series titles:", paste(titles, collapse = " | "), "\n")

# match file prefixes (GSMxxxxxxxx_...) to titles
pre_gsm <- sub("_.*$", "", prefixes)
cond <- titles[match(pre_gsm, gsms)]
if (any(is.na(cond))) { cat("WARNING: GSM match failed, falling back to prefix parsing\n") }
print(data.frame(prefix = prefixes, condition = cond))

parse_cond <- function(tt) {
  # "Liver, Peak CCl4 Fibrosis" / "Liver, 2wk CCl4 Regression" / "Liver, Peak FAT-MASH Regression" ...
  model <- ifelse(grepl("FAT", tt), "FAT-MASH", "CCl4")
  tp <- dplyr::case_when(grepl("Peak", tt) ~ "peak",
                         grepl("2wk", tt) ~ "reg2wk",
                         grepl("6wk", tt) ~ "reg6wk", TRUE ~ "other")
  list(model = model, tp = tp)
}
pc <- lapply(cond, parse_cond)

objs <- lapply(seq_along(mats), function(i) {
  o <- CreateSeuratObject(mats[[i]], project = prefixes[i], min.cells = 3, min.features = 200)
  o$sample <- prefixes[i]; o$model <- pc[[i]]$model; o$timepoint <- pc[[i]]$tp
  o$condition <- paste(pc[[i]]$model, pc[[i]]$tp, sep = "_")
  o
})

## ---- QC ----
qc_stats <- lapply(objs, function(o) {
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^mt-")
  data.frame(sample = o$sample[1], n_cells_raw = ncol(o),
             median_nFeature = median(o$nFeature_RNA), median_nCount = median(o$nCount_RNA),
             median_mt = median(o$percent.mt))
}) |> bind_rows()
print(qc_stats)
write.csv(qc_stats, file.path(OUT, "qc_mouse_sc_raw.csv"), row.names = FALSE)

objs <- lapply(objs, function(o) {
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^mt-")
  subset(o, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 15)
})
cat("cells after QC:\n"); print(sapply(objs, ncol))

## ---- doublets (scDblFinder) ----
suppressMessages(library(scDblFinder))
objs <- lapply(objs, function(o) {
  sce <- as.SingleCellExperiment(o)
  sce <- scDblFinder(sce, dbr.sd = 1)
  o$doublet <- sce$scDblFinder.class
  o[, o$doublet == "singlet"]
})
cat("cells after doublet removal:\n"); print(sapply(objs, ncol))

## ---- integrate ----
objs <- lapply(objs, function(o) {
  o <- NormalizeData(o, verbose = FALSE)
  o <- FindVariableFeatures(o, nfeatures = 3000, verbose = FALSE)
  o
})
feats_int <- SelectIntegrationFeatures(objs, nfeatures = 3000)
m <- merge(objs[[1]], objs[-1])
VariableFeatures(m) <- feats_int
m <- ScaleData(m, verbose = FALSE)
m <- RunPCA(m, npcs = 50, verbose = FALSE)
suppressMessages(library(harmony))
m <- RunHarmony(m, group.by.vars = "sample", verbose = FALSE)
m <- RunUMAP(m, reduction = "harmony", dims = 1:30, verbose = FALSE)
m <- FindNeighbors(m, reduction = "harmony", dims = 1:30, verbose = FALSE)
m <- FindClusters(m, resolution = 0.5, verbose = FALSE)
saveRDS(m, file.path(OUT, "mouse_sc_integrated.rds"))

## ---- annotate by canonical markers ----
markers <- list(
  Hepatocyte = c("Alb", "Apoa2", "Ttr", "Serpina1c"),
  HSC_myofibroblast = c("Col1a1", "Col3a1", "Acta2", "Pdgfrb", "Des", "Rgs5"),
  Endothelial_LSEC = c("Pecam1", "Kdr", "Stab2", "Clec4g", "Flt4"),
  Kupffer = c("Marco", "Cd5l", "Vsig4", "Timd4", "Clec4f"),
  SAMac_Trem2 = c("Trem2", "Cd9", "Spp1", "Gpnmb", "Fabp5", "Lgals3"),
  Monocyte_MF = c("Ly6c2", "Plac8", "Chil3", "Ccr2"),
  T_NK = c("Cd3d", "Cd3e", "Cd8a", "Nkg7"),
  B_cell = c("Cd79a", "Ms4a1", "Cd19"),
  Cholangiocyte = c("Krt19", "Krt7", "Epcam", "Sox9"),
  Dendritic = c("Itgax", "Clec9a", "Batf3", "Siglech")
)
suppressMessages(library(UCell))
m <- AddModuleScore_UCell(m, features = markers, name = "_mark")
score_cols <- paste0(names(markers), "_mark")
stopifnot(all(score_cols %in% colnames(m@meta.data)))
cl_avg <- sapply(score_cols, function(v) tapply(m@meta.data[[v]], m$seurat_clusters, mean, na.rm = TRUE))
# cl_avg: rows = clusters, cols = marker sets
ct_map <- setNames(colnames(cl_avg)[apply(cl_avg, 1, which.max)], rownames(cl_avg))
m$celltype <- ct_map[as.character(m$seurat_clusters)]
print(table(m$celltype, m$condition))
write.csv(as.data.frame(cl_avg), file.path(OUT, "mouse_cluster_marker_scores.csv"))
write.csv(m@meta.data, file.path(OUT, "mouse_sc_metadata.csv"))
saveRDS(ct_map, file.path(OUT, "mouse_celltype_map.rds"))
cat("DONE\n")
