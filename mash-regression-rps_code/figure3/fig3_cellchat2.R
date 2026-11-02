#!/usr/bin/env Rscript
# CellChat part 2: manual differential analysis from saved per-condition objects
suppressMessages({ library(CellChat); library(dplyr); library(tidyr) })
OUT <- "/workspace/fig3_out"
ccs <- readRDS(file.path(OUT, "cellchat_objects.rds"))

melt_net <- function(mat, measure, cond, model) {
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("source", "target", "value")
  df$measure <- measure; df$condition <- cond; df$model <- model
  df
}

all_pairs <- list(); all_path <- list()
for (mdl in c("CCl4", "FAT-MASH")) {
  for (tp in c("peak", "reg6wk")) {
    cc <- ccs[[paste0(mdl, "_", tp)]]
    all_pairs[[paste(mdl, tp)]] <- melt_net(cc@net$count, "count", tp, mdl)
    all_pairs[[paste(mdl, tp, "w")]] <- melt_net(cc@net$weight, "weight", tp, mdl)
    # pathway-level total strength
    pw <- apply(cc@netP$prob, 3, sum, na.rm = TRUE)
    all_path[[paste(mdl, tp)]] <- data.frame(model = mdl, condition = tp,
                                             pathway = names(pw), strength = pw)
  }
}
pairs <- bind_rows(all_pairs)
paths <- bind_rows(all_path)

# differential: reg6wk vs peak per model
diff_pairs <- pairs %>%
  pivot_wider(names_from = condition, values_from = value, values_fill = 0) %>%
  mutate(delta = reg6wk - peak)
write.csv(diff_pairs, file.path(OUT, "cellchat_pair_deltas.csv"), row.names = FALSE)

diff_paths <- paths %>%
  pivot_wider(names_from = condition, values_from = strength, values_fill = 0) %>%
  mutate(delta = reg6wk - peak, log2fc = log2((reg6wk + 1e-6) / (peak + 1e-6)))
write.csv(diff_paths, file.path(OUT, "cellchat_pathway_deltas.csv"), row.names = FALSE)

cat("== pathways up in regression (both models) ==\n")
up <- diff_paths %>% group_by(pathway) %>% summarise(n_up = sum(delta > 0), mean_delta = mean(delta), .groups = "drop") %>%
  filter(n_up == 2) %>% arrange(desc(mean_delta))
print(head(up, 20))
cat("== pathways down in regression (both models) ==\n")
dn <- diff_paths %>% group_by(pathway) %>% summarise(n_dn = sum(delta < 0), mean_delta = mean(delta), .groups = "drop") %>%
  filter(n_dn == 2) %>% arrange(mean_delta)
print(head(dn, 20))

# myeloid -> HSC LR tables per condition
mye <- c("SAMac (Trem2+)", "Monocyte/MF", "Kupffer", "Neutrophil", "DC")
for (mdl in c("CCl4", "FAT-MASH")) {
  for (tp in c("peak", "reg6wk")) {
    cc <- ccs[[paste0(mdl, "_", tp)]]
    df <- subsetCommunication(cc)
    df2 <- df[df$source %in% mye & df$target == "HSC/myofibroblast", ]
    write.csv(df2, file.path(OUT, paste0("cellchat_myeloid_to_hsc_", mdl, "_", tp, ".csv")), row.names = FALSE)
    cat(mdl, tp, "myeloid->HSC interactions:", nrow(df2), "\n")
  }
}
cat("DONE\n")
