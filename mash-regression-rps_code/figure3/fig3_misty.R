#!/usr/bin/env Rscript
# Figure 3 panel f (mouse part): MISTy — which cell-type neighborhoods predict
# the reversal/emergent spot scores on Visium (per slide)
suppressMessages({ library(mistyR); library(dplyr); library(tidyr) })
OUT <- "/workspace/fig3_out"
DATA <- "/workspace/fig3_data/GSE306327"

read_pos <- function(posdir) {
  pos <- read.csv(file.path(DATA, posdir, "spatial", "tissue_positions.csv"))
  data.frame(barcode = pos$barcode, x = pos$pxl_col_in_fullres, y = pos$pxl_row_in_fullres)
}

for (sn in c("peak", "regression")) {
  w <- read.csv(file.path(OUT, paste0("rctd_weights_", sn, ".csv")), row.names = 1, check.names = FALSE)
  sc <- read.csv(file.path(OUT, paste0("spot_scores_", sn, ".csv")))
  pos <- read_pos(if (sn == "peak") "peak" else "regression")
  df <- pos %>% inner_join(sc, by = "barcode") %>%
    inner_join(w %>% tibble::rownames_to_column("barcode"), by = "barcode")
  df <- df[!is.na(df$reversal_prog), ]
  cat(sn, "spots with all data:", nrow(df), "\n")

  expr <- df %>% select(reversal_prog, emergent_prog, all_of(colnames(w))) %>% as.data.frame()
  colnames(expr) <- make.names(colnames(expr))  # ranger-safe names
  coords <- df %>% select(x, y)
  # standardize
  expr <- as.data.frame(scale(expr))

  view <- create_initial_view(expr) %>% add_paraview(coords, l = 200)  # ~2 spot diameters
  run_misty(view, results.folder = file.path(OUT, paste0("misty_", sn)), cached = FALSE)
  res <- collect_results(file.path(OUT, paste0("misty_", sn)))
  imp <- res$importances.aggregated
  write.csv(imp, file.path(OUT, paste0("misty_importances_", sn, ".csv")), row.names = FALSE)
  # R2
  write.csv(res$improvements, file.path(OUT, paste0("misty_improvements_", sn, ".csv")), row.names = FALSE)
  cat(sn, "done\n")
}
cat("DONE\n")
