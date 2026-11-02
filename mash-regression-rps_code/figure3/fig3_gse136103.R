#!/usr/bin/env Rscript
# Figure 3 supplement: GSE136103 (Ramachandran 2019 human cirrhosis atlas) cross-validation
# of component localization. cd45+ and cd45- fractions per donor; healthy vs cirrhotic.
suppressMessages({ library(Seurat); library(dplyr); library(Matrix); library(UCell) })
DATA <- "/workspace/fig3_data/GSE136103"; OUT <- "/workspace/fig3_out"

feats <- list.files(DATA, pattern = "_genes.tsv.gz$")
prefixes <- sub("_genes.tsv.gz$", "", feats)
prefixes <- prefixes[grepl("healthy|cirrhotic", prefixes)]   # liver only; drop blood + mouse
cat("samples:", length(prefixes), "\n")

read_sample <- function(px) {
  feat <- read.delim(gzfile(file.path(DATA, paste0(px, "_genes.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)
  mat <- readMM(gzfile(file.path(DATA, paste0(px, "_matrix.mtx.gz"))))
  bc <- read.delim(gzfile(file.path(DATA, paste0(px, "_barcodes.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)[[1]]
  colnames(mat) <- bc
  rownames(mat) <- make.unique(if (ncol(feat) >= 2) feat[[2]] else feat[[1]])
  t(mat)  # older 10X: genes x cells already? check dims: matrix.mtx header rows=genes
}
# verify orientation on first sample
m0 <- readMM(gzfile(file.path(DATA, paste0(prefixes[1], "_matrix.mtx.gz"))))
cat("first matrix dims:", nrow(m0), "x", ncol(m0), "\n")

mats <- lapply(prefixes, function(px) {
  feat <- read.delim(gzfile(file.path(DATA, paste0(px, "_genes.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)
  mat <- readMM(gzfile(file.path(DATA, paste0(px, "_matrix.mtx.gz"))))
  bc <- read.delim(gzfile(file.path(DATA, paste0(px, "_barcodes.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)[[1]]
  if (nrow(mat) == length(bc)) mat <- t(mat)   # cells x genes -> genes x cells
  colnames(mat) <- bc
  rownames(mat) <- make.unique(if (ncol(feat) >= 2) feat[[2]] else feat[[1]])
  mat
})
names(mats) <- prefixes

meta <- data.frame(sample = prefixes) %>%
  mutate(donor = sub("^GSM[0-9]+_((healthy|cirrhotic)[0-9]+)_.*$", "\\1", sample),
         group = ifelse(grepl("healthy", sample), "Healthy", "Cirrhotic"),
         fraction = ifelse(grepl("cd45\\+", sample), "CD45pos", "CD45neg"))
print(meta)

objs <- lapply(seq_along(mats), function(i) {
  o <- CreateSeuratObject(mats[[i]], project = prefixes[i], min.cells = 3, min.features = 200)
  o$sample <- prefixes[i]; o$donor <- meta$donor[i]; o$group <- meta$group[i]; o$fraction <- meta$fraction[i]
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  subset(o, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 15)
})
nc <- sapply(objs, ncol); print(nc)
objs <- objs[nc >= 50]; cat("samples kept:", length(objs), "\n")
m <- merge(objs[[1]], objs[-1])
m <- JoinLayers(m)
m <- NormalizeData(m, verbose = FALSE)
m <- FindVariableFeatures(m, nfeatures = 2000, verbose = FALSE)
m <- ScaleData(m, verbose = FALSE)
m <- RunPCA(m, npcs = 30, verbose = FALSE)
nd <- min(20, ncol(Embeddings(m, "pca")))
cat("PCA dims available:", ncol(Embeddings(m, "pca")), " using:", nd, "\n")
suppressMessages(library(harmony))
m <- RunHarmony(m, group.by.vars = "sample", dims.use = 1:nd, verbose = FALSE)
m <- RunUMAP(m, reduction = "harmony", dims = 1:nd, verbose = FALSE)
m <- FindNeighbors(m, reduction = "harmony", dims = 1:nd, verbose = FALSE) |> FindClusters(resolution = 0.4, verbose = FALSE)

markers <- list(
  Hepatocyte = c("ALB", "APOA2", "TTR"),
  HSC_myofibroblast = c("COL1A1", "COL3A1", "ACTA2", "PDGFRB"),
  Endothelial_LSEC = c("PECAM1", "STAB2", "CLEC4G"),
  Kupffer = c("MARCO", "CD5L", "VSIG4", "TIMD4"),
  SAMac_OLR1 = c("OLR1", "TREM2", "CD9", "SPP1", "GPNMB"),
  Monocyte_MF = c("LYZ", "S100A8", "S100A9", "CCR2"),
  T_NK = c("CD3D", "CD3E", "NKG7"),
  B_cell = c("CD79A", "MS4A1"),
  Plasma = c("MZB1", "JCHAIN"),
  Cholangiocyte = c("KRT19", "KRT7", "EPCAM"),
  Dendritic = c("ITGAX", "CLEC9A")
)
d <- GetAssayData(m, layer = "data")
sc <- as.data.frame(ScoreSignatures_UCell(d, features = markers, name = "_mark"))
stopifnot(identical(rownames(sc), colnames(m)))
m@meta.data <- cbind(m@meta.data, sc)
score_cols <- colnames(sc)
cl_avg <- sapply(score_cols, function(v) tapply(m@meta.data[[v]], m$seurat_clusters, mean, na.rm = TRUE))
zs <- scale(cl_avg)
ct_map <- setNames(colnames(zs)[apply(zs, 1, which.max)], rownames(zs))
m$celltype <- unname(ct_map[as.character(m$seurat_clusters)])
print(table(m$celltype, m$group))

# component scoring
rev_g <- intersect(readLines(file.path(OUT, "genes_reversal_human.txt")), rownames(m))
eme_g <- intersect(readLines(file.path(OUT, "genes_emergent_human.txt")), rownames(m))
scp <- as.data.frame(ScoreSignatures_UCell(d, features = list(reversal = rev_g, emergent = eme_g), name = "_prog"))
m@meta.data <- cbind(m@meta.data, scp)
saveRDS(m, file.path(OUT, "gse136103_scored_tmp.rds"))
file.rename(file.path(OUT, "gse136103_scored_tmp.rds"), file.path(OUT, "gse136103_scored.rds"))
write.csv(m@meta.data, file.path(OUT, "gse136103_metadata.csv"))

tests <- bind_rows(lapply(c("reversal", "emergent"), function(comp) {
  v <- paste0(comp, "_prog")
  bind_rows(lapply(unique(m$celltype), function(ct) {
    sub <- m@meta.data[m@meta.data$celltype == ct, ]
    wt <- tryCatch(wilcox.test(sub[[v]] ~ sub$group)$p.value, error = function(e) NA_real_)
    data.frame(comp = comp, celltype = ct, p = wt,
               med_healthy = median(sub[[v]][sub$group == "Healthy"]),
               med_cirrhotic = median(sub[[v]][sub$group == "Cirrhotic"]))
  }))
}))
tests$fdr <- p.adjust(tests$p, "BH")
write.csv(tests, file.path(OUT, "gse136103_ucell_tests.csv"), row.names = FALSE)
print(tests %>% arrange(fdr) %>% head(10))
cat("DONE\n")
