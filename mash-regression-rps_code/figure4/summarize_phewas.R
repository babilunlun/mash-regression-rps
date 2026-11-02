.libPaths(c("/workspace/fig4/Rlib", .libPaths()))
suppressMessages(library(data.table))
pw <- fread("/workspace/fig4/qc/phewas_associations.csv")
pw[, p := as.numeric(p)]
cat("associations:", nrow(pw), "| SNPs with hits:", length(unique(pw$rsid)),
    "| genes:", length(unique(pw$gene)), "\n")

# categorize: self-pQTL vs clinical
pw[, cat := fcase(
  grepl("measurement|protein|level|ratio|amount", trait, ignore.case = TRUE), "molecular/protein",
  grepl("count|percentage|volume|hemoglobin|platelet|leukocyte|lymphocyte|neutrophil|monocyte|eosinophil|basophil|erythrocyte|reticulocyte|blood", trait, ignore.case = TRUE), "blood trait",
  grepl("cholesterol|lipid|triglycer|lipoprotein|apolipo", trait, ignore.case = TRUE), "lipid",
  grepl("diabet|glucose|insulin|body mass|obesity|weight|waist", trait, ignore.case = TRUE), "metabolic",
  default = "clinical/other")]
print(pw[, .N, by = cat][order(-N)])

cat("\n=== clinical/other associations (safety-relevant) ===\n")
clin <- pw[cat == "clinical/other"][order(p)]
print(clin[, .(gene, rsid, trait, p = signif(p, 2))][1:min(25, nrow(clin))])

cat("\n=== per-gene strongest non-molecular association ===\n")
nm <- pw[cat != "molecular/protein", .(minp = min(p, na.rm = TRUE),
                                       trait = trait[which.min(p)]), by = .(gene, rsid)]
print(nm[order(minp)][1:min(20, nrow(nm))])
fwrite(pw, "/workspace/fig4/qc/phewas_associations_categorized.csv", quote = TRUE)
