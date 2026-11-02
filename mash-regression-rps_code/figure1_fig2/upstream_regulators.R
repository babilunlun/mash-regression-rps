#!/usr/bin/env Rscript
## upstream_regulators.R v2 — TF-target enrichment (ChEA_2022 + TRRUST_2019, local GMT)
## TF-level sets (union across experiments), universe = 27,265 genes tested in paired DE
suppressMessages({library(data.table); library(clusterProfiler); library(ggplot2)})
OUT <- "/workspace/upstream"; dir.create(OUT, showWarnings = FALSE)

emg <- trimws(readLines("/mnt/shared-workspace/shared/human_rebuild/genes_emergent_human.txt")); emg <- emg[nzchar(emg)]
rev <- trimws(readLines("/mnt/shared-workspace/shared/human_rebuild/genes_reversal_human.txt")); rev <- rev[nzchar(rev)]
de <- fread("/mnt/results/Figure1/tables/de_paired41.csv")
universe <- de$gene
down250 <- de[adj.P.Val < 0.05 & logFC < 0, gene]
cat("emergent:", length(emg), "| reversal:", length(rev), "| down250:", length(down250), "| universe:", length(universe), "\n")

read_gmt_tf <- function(path, source) {
  ln <- readLines(path)
  rbindlist(lapply(ln, function(l) {
    f <- strsplit(l, "\t")[[1]]
    tf <- toupper(strsplit(f[1], "[ _]")[[1]][1])
    genes <- f[-(1:2)]; genes <- genes[nzchar(genes)]
    data.table(tf = paste0(source, ":", tf), gene_symbol = genes)
  }))
}
tft <- unique(rbind(read_gmt_tf("/workspace/chea2022.gmt", "ChEA2022"),
                    read_gmt_tf("/workspace/trrust2019.gmt", "TRRUST2019")))
tft <- tft[gene_symbol %in% universe]
cat("TF sets:", length(unique(tft$tf)), "| T2G pairs:", nrow(tft), "\n")

run_enr <- function(genes, name) {
  e <- enricher(genes, universe = universe, TERM2GENE = as.data.frame(tft), pAdjustMethod = "BH",
                pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
  d <- as.data.table(e)
  if (nrow(d) == 0) return(d)
  d[, list_name := name]
  d
}
res <- rbind(run_enr(emg, "emergent_214"), run_enr(rev, "reversal_65"), run_enr(down250, "down_250"))
fwrite(res, file.path(OUT, "tf_enrichment_chea_trrust.csv"))

for (nm in c("emergent_214", "reversal_65", "down_250")) {
  d <- res[list_name == nm][order(p.adjust)]
  cat("\n== ", nm, " ==  FDR<0.05:", sum(d$p.adjust < 0.05), "of", nrow(d), "\n")
  print(d[1:min(12, .N), .(ID, Count, pvalue, p.adjust)])
  atf4 <- d[grepl("ATF4", ID)]
  if (nrow(atf4)) cat("ATF4: ", paste(sprintf("%s FDR=%.3g", atf4$ID, atf4$p.adjust), collapse = "; "), "\n")
}

## dotplot: top 10 per list
top <- res[, .SD[order(p.adjust)][1:min(10, .N)], by = list_name]
or_tab <- top[, {
  parts <- do.call(rbind, strsplit(GeneRatio, "/")); bg <- do.call(rbind, strsplit(BgRatio, "/"))
  .(OR = (as.numeric(parts[,1]) / as.numeric(parts[,2])) / (as.numeric(bg[,1]) / as.numeric(bg[,2])))
}, by = .(list_name, ID)]
top <- merge(top, or_tab, by = c("list_name", "ID"))
top[, list_name := factor(list_name, levels = c("emergent_214", "reversal_65", "down_250"),
                          labels = c("Emergent (214)", "Reversal (65)", "Down in regression (250)"))]
tf_order <- top[, .(minfdr = min(p.adjust)), by = ID][order(minfdr)]$ID
top[, ID := factor(ID, levels = rev(unique(tf_order)))]
p <- ggplot(top, aes(list_name, ID)) +
  geom_point(aes(size = OR, color = p.adjust)) +
  scale_color_gradient(low = "#D55E00", high = "#0072B2", trans = "log10", name = "FDR",
                       guide = guide_colorbar(barwidth = 0.6, barheight = 4)) +
  scale_size_continuous(name = "Odds ratio", range = c(0.8, 5)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9, base_family = "Liberation Sans") +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1),
        axis.text.y = element_text(size = 7))
ggsave(file.path(OUT, "Fig2-S9_upstream_regulators.svg"), p, width = 4.8, height = 5.6, device = svglite::svglite)
ggsave(file.path(OUT, "Fig2-S9_upstream_regulators.png"), p, width = 4.8, height = 5.6, dpi = 400)
cat("UPSTREAM_DONE\n")
