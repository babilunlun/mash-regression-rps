#!/usr/bin/env Rscript
# install_r_pkgs.R — minimal R package set for Figure 4
.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
dir.create("/workspace/fig4/Rlib", showWarnings = FALSE)
pkgs <- c("data.table", "ggplot2", "svglite", "patchwork", "coloc", "scales")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE, lib.loc = "/workspace/fig4/Rlib")) {
    message("installing ", p)
    install.packages(p, lib = "/workspace/fig4/Rlib", quiet = TRUE, Ncpus = 4)
  }
}
# coloc wants susieR optionally; skip. Verify:
for (p in pkgs) cat(p, requireNamespace(p, quietly = TRUE, lib.loc = "/workspace/fig4/Rlib"), "\n")
cat("R_PKGS_DONE\n")
