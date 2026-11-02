# MASH fibrosis regression program and the Regression Program Score (RPS)

Analysis code accompanying the manuscript:

**A conserved transcriptional program of MASH fibrosis regression yields a sample-level score for treatment response** (under review).

## Contents

| Folder | Contents |
|---|---|
| `figure1_fig2/` | Bulk RNA-seq differential expression, replication, intervention-stratified analyses, upstream-regulator enrichment (Figs. 1–2) |
| `figure3/` | Single-cell (human/mouse), spatial (Visium RCTD, GeoMx, CosMx, MISTy), CellChat, SAMac analyses, program scoring (Fig. 3) |
| `figure4/` | GWAS download/build, Mendelian randomization, coloc/SuSiE, PheWAS, evidence-tier meta-analysis (Fig. 4) |
| `figure5/` | Target prioritization: Open Targets tractability, gnomAD constraint, GTEx specificity, TP53INP1 analyses (Fig. 5) |
| `figure6/` | RPS construction, longitudinal validation, sensitivity and leave-one-out analyses (Fig. 6) |
| `figure7/` | Clinical benchmarking: diagnostic groups, BARD comparison, baseline-F3 stratum statistics, PLS-NAFLD bridge, secreted subset (Fig. 7, Fig. S7) |
| `rps_calculator/` | Source code of the interactive RPS calculator (Shiny) |
| `sessionInfo/` | R session information for the main analysis phase and the v2 revision phase |

See `MANIFEST.md` for a per-script table (in Chinese).

## Notes

- R 4.4.2; a small number of helper scripts run under python3.
- Scripts are organized by the main figure they support; `*_v2` / `*_v3` suffixes denote iterative versions. The final output version of each figure is recorded in the corresponding figure README of the manuscript repository.
- Paths inside scripts refer to the original analysis workspace (`/workspace/...`) and must be adapted to your local layout. Input data are public: GEO (GSE106737, GSE83452, GSE135251, GSE136103, GSE89632, GSE193066, GSE269493, GSE294630, GSE70779, GSE298719, GSE305330, GSE287826, GSE292268, GSE289173), eQTLGen, GTEx, INTERVAL, FinnGen, GWAS Catalog (GCST90091033, GCST90267352, GCST90590317), gnomAD, Open Targets.
- The RPS is computed from within-sample gene ranks as the difference of two gene-set activities (up arm minus down arm); it requires neither cohort-level normalization nor an outcome-trained model and can be applied to a single new sample on any transcriptomic platform. See `figure6/fig6_prep_score.R` and `figure6/fig6_score_gse83452.R`.
- Three steps are archived as intermediate data or final artwork rather than as scripts: the builders of the two Fig. 7 input tables and of the secretome-annotation table (the tables themselves are included in `figure7/`), and the plotting scripts for the Fig. 7 PLS-bridge and secretome panels (final artwork is deposited with the manuscript figures). The calculator's reference matrix (`reference_expr.rds`, 7.6 MB binary) and example input file are not included; they can be rebuilt from public data (GSE193066) with the `figure6/` scoring scripts.

## License

Code is provided for academic research use.
