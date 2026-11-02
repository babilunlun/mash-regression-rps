#!/usr/bin/env Rscript
## fig1_stratified.R — A5: RYGB vs lifestyle stratified DE; A6: responder vs non-responder baseline
suppressMessages({library(data.table); library(limma); library(GSVA)})
OUT <- "/workspace/fig1_strat"; dir.create(OUT, showWarnings = FALSE)

## ---------- expression from series matrix (RMA+SVA by submitter) ----------
f <- "/workspace/geo/GSE106737_series_matrix.txt.gz"
ln <- readLines(f, n = 300); beg <- grep("!series_matrix_table_begin", ln)
d <- fread(f, skip = beg, header = TRUE)
suppressMessages(library(hugene20sttranscriptcluster.db))
map <- AnnotationDbi::select(hugene20sttranscriptcluster.db, keys = as.character(d$ID_REF),
                             keytype = "PROBEID", columns = "SYMBOL")
map <- as.data.table(map)[!is.na(SYMBOL)]; map <- map[!duplicated(PROBEID)]
d2 <- d[as.character(ID_REF) %in% map$PROBEID]
sym <- map$SYMBOL[match(as.character(d2$ID_REF), map$PROBEID)]
m <- as.matrix(d2[, -1]); if (max(m) > 100) m <- log2(m + 1)
o <- order(rowMeans(m), decreasing = TRUE); keep <- o[!duplicated(sym[o])]
m <- m[keep, ]; rownames(m) <- sym[keep]
cat("gene-level matrix:", nrow(m), "x", ncol(m), "\n")

## ---------- metadata ----------
pg <- read.table("/workspace/geo/gse106737_pg.txt", col.names = c("gsm", "pgroup"))
meta <- as.data.table(merge(read.csv("/mnt/results/Figure1/tables/metadata_GSE106737.csv"), pg, by = "gsm"))
bd <- read.csv("/mnt/results/Figure1/tables/blue_delta_by_pair.csv")
meta[, arm := fcase(pgroup == 1, "RYGB", pgroup %in% c(2, 3), "lifestyle", default = NA_character_)]
paired <- meta[sample_type == "paired"]
paired <- merge(paired, bd[, c("pair_id", "reconstructed_response")], by = "pair_id")
cat("arms:", table(unique(paired[, .(pair_id, arm)])$arm), "\n")

de41 <- fread("/mnt/results/Figure1/tables/de_paired41.csv")
prog353 <- de41[adj.P.Val < 0.05]$gene
cat("program genes (FDR<0.05):", length(prog353), "\n")

run_paired_de <- function(gsms, pair_ids) {
  mm <- m[, gsms]
  ph <- data.frame(row.names = gsms, time = factor(ifelse(grepl("baseline", paired$title[match(gsms, paired$gsm)]), "base", "fu")),
                   pair = factor(pair_ids))
  des <- model.matrix(~ time + pair, ph)
  fit <- eBayes(lmFit(mm, des))
  tt <- topTable(fit, coef = "timefu", number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, t = tt$t, P = tt$P.Value, FDR = tt$adj.P.Val)
}

## ---------- sanity: all-pairs DE vs archived de_paired41 ----------
all_de <- run_paired_de(paired$gsm, paired$pair_id)
chk <- merge(all_de, de41, by = "gene", suffixes = c("", "_arch"))
cat("sanity all-pairs vs archived: Spearman", round(cor(chk$logFC, chk$logFC_arch, method = "spearman"), 4),
    "| sign concordance in 353:", with(chk[gene %in% prog353],
    sum(sign(logFC) == sign(logFC_arch))), "/", length(prog353), "\n")

## ---------- A5: arm-stratified DE ----------
ry <- paired[arm == "RYGB"]; ls_ <- paired[arm == "lifestyle"]
de_rygb <- run_paired_de(ry$gsm, ry$pair_id)
de_life <- run_paired_de(ls_$gsm, ls_$pair_id)
fwrite(de_rygb, file.path(OUT, "de_stratified_rygb.csv"))
fwrite(de_life, file.path(OUT, "de_stratified_lifestyle.csv"))

cmp <- merge(de_rygb, de_life, by = "gene", suffixes = c("_rygb", "_life"))
cmp353 <- cmp[gene %in% prog353]
rho353 <- cor(cmp353$logFC_rygb, cmp353$logFC_life, method = "spearman")
sign353 <- sum(sign(cmp353$logFC_rygb) == sign(cmp353$logFC_life))
rho_all <- cor(cmp$logFC_rygb, cmp$logFC_life, method = "spearman")
cat(sprintf("A5: 353 genes RYGB-vs-lifestyle rho=%.3f, sign %d/353; genome-wide rho=%.3f\n",
            rho353, sign353, rho_all))
## per-arm program replication: how many of 353 are nominally significant & same direction in each arm
de41s <- de41[match(prog353, gene)]
for (nm in c("rygb", "life")) {
  dd <- get(paste0("de_", nm))
  mg <- merge(de41s[, .(gene, prog_logFC = logFC)], dd, by = "gene")
  cat(sprintf("  %s: nominal P<0.05 & same direction: %d/353; median |logFC| %.3f\n", nm,
              sum(mg$P < 0.05 & sign(mg$logFC) == sign(mg$prog_logFC)), median(abs(mg$logFC))))
}
## interaction via per-pair deltas (paired design: delta ~ arm IS the time x arm interaction)
bp <- paired[timepoint == "baseline", .(pair_id, gsm_b = gsm, arm, sex)]
fp <- paired[timepoint == "followup", .(pair_id, gsm_f = gsm)]
dmm <- merge(bp, fp, by = "pair_id")
dmm$arm <- factor(dmm$arm, levels = c("RYGB", "lifestyle"))
delta <- m[, dmm$gsm_f] - m[, dmm$gsm_b]
colnames(delta) <- dmm$pair_id
des_d <- model.matrix(~ arm + sex, data = dmm)
fit_d <- eBayes(lmFit(delta, des_d))
tt_i <- topTable(fit_d, coef = "armlifestyle", number = Inf, sort.by = "none")
int <- data.table(gene = rownames(tt_i), logFC_int = tt_i$logFC, P_int = tt_i$P.Value, FDR_int = tt_i$adj.P.Val)
fwrite(int, file.path(OUT, "de_interaction_arm.csv"))
cat("A5 interaction: FDR<0.05 genome-wide:", sum(int$FDR_int < 0.05),
    "| among 353:", sum(int$FDR_int < 0.05 & int$gene %in% prog353),
    "| nominal among 353:", sum(int$P_int < 0.05 & int$gene %in% prog353), "\n")

## ---------- A6: responder vs non-responder baseline ----------
base <- paired[timepoint == "baseline"]
cat("A6 baseline: R", sum(base$reconstructed_response == "responder"),
    "NR", sum(base$reconstructed_response == "nonresponder"), "\n")
phb <- data.frame(row.names = base$gsm,
                  resp = factor(base$reconstructed_response, levels = c("nonresponder", "responder")),
                  sex = factor(base$sex), arm = factor(base$arm))
desb <- model.matrix(~ resp + sex + arm, phb)  # all paired patients are NASH+T2D: t2d not estimable
fitb <- eBayes(lmFit(m[, base$gsm], desb))
tt_b <- topTable(fitb, coef = "respresponder", number = Inf, sort.by = "none")
bde <- data.table(gene = rownames(tt_b), logFC = tt_b$logFC, t = tt_b$t, P = tt_b$P.Value, FDR = tt_b$adj.P.Val)
fwrite(bde, file.path(OUT, "baseline_R_vs_NR_de.csv"))
cat("A6: baseline R-vs-NR FDR<0.05 genes:", sum(bde$FDR < 0.05), "| nominal:", sum(bde$P < 0.05), "\n")
## do program genes differ at baseline? enrichment of nominal P among 353
bd353 <- bde[gene %in% prog353]
cat("A6: 353 genes nominal-P<0.05 at baseline:", sum(bd353$P < 0.05),
    sprintf("(%.1f%% vs genome %.1f%%)", 100*mean(bd353$P < 0.05), 100*mean(bde$P < 0.05)), "\n")
print(wilcox.test(bd353$P, bde$P))
## baseline RPS
sg <- fread("/mnt/results/Figure6/tables/score_genes.csv")
par <- ssgseaParam(exprData = m[, base$gsm],
                   geneSets = list(up = intersect(sg[arm == "up"]$gene, rownames(m)),
                                   down = intersect(sg[arm == "down"]$gene, rownames(m))))
sc <- gsva(par, verbose = FALSE)
rps <- data.table(gsm = colnames(m[, base$gsm]), RPS = sc["up", ] - sc["down", ])
rps <- merge(rps, base[, .(gsm, reconstructed_response, arm)], by = "gsm")
w <- wilcox.test(RPS ~ reconstructed_response, rps)
aucv <- { r <- rank(rps$RPS); y <- rps$reconstructed_response == "responder"
          (sum(r[y]) - sum(y)*(sum(y)+1)/2) / (sum(y)*sum(!y)) }
cat(sprintf("A6: baseline RPS R vs NR: Wilcoxon P = %.3f, AUC = %.3f\n", w$p.value, aucv))
cat("   median RPS by group:", rps[, median(RPS), by = reconstructed_response]$V1, "\n")
fwrite(rps, file.path(OUT, "baseline_rps_response.csv"))

summ <- data.table(
  metric = c("sanity_rho_allpairs_vs_archived", "sign_concordance_353_vs_archived",
             "rho353_rygb_vs_life", "sign353_rygb_vs_life", "rho_genomewide_rygb_vs_life",
             "interaction_FDR05_genomewide", "interaction_FDR05_in353", "interaction_nominal_in353",
             "baseline_RNR_FDR05", "baseline_RNR_nominal", "baseline_353_nominal_frac",
             "baseline_RPS_wilcox_P", "baseline_RPS_AUC"),
  value = c(round(cor(chk$logFC, chk$logFC_arch, method = "spearman"), 4),
            paste0(with(chk[gene %in% prog353], sum(sign(logFC) == sign(logFC_arch))), "/353"),
            round(rho353, 3), paste0(sign353, "/353"), round(rho_all, 3),
            sum(int$FDR_int < 0.05), sum(int$FDR_int < 0.05 & int$gene %in% prog353),
            sum(int$P_int < 0.05 & int$gene %in% prog353),
            sum(bde$FDR < 0.05), sum(bde$P < 0.05), round(mean(bd353$P < 0.05), 3),
            signif(w$p.value, 3), round(aucv, 3)))
fwrite(summ, file.path(OUT, "stratified_summary.csv"))
print(summ)
cat("STRAT_DONE\n")
