# RPS calculator — Regression Potential Score online tool

A Shiny app that computes the **Regression Potential Score (RPS)** for a single
liver expression sample and places it in reference cohorts of MASH/NAFLD
biopsies.

**RPS = ssGSEA(up arm, 22 genes) − ssGSEA(down arm, 138 genes)**

The two arms are the strict cross-species conserved core of the MASH fibrosis
regression program (160 genes). Higher RPS indicates an expression state closer
to fibrosis regression; within-patient RPS increases track histological
improvement.

## Input

A CSV file with two columns: **gene symbol** and **normalized expression**
(e.g. logCPM, log2-normalized counts, or microarray intensity). One sample per
upload. Because ssGSEA uses within-sample ranks, any monotone transformation of
the values yields the identical score.

`example_input.csv` is a real baseline biopsy profile (GSE193066 sample
HUnafld001) that can be loaded from the app to try it out; it reproduces the
published RPS of that sample exactly (0.6324901).

## Output

- RPS and its up/down arm components
- Percentile against a selectable reference panel:
  - GSE193066 baseline (n = 106) and follow-up (n = 58), bariatric cohort, RNA-seq
  - GSE135251 NAFLD spectrum (n = 216), RNA-seq
  - GSE89632 NAFLD spectrum (n = 56), microarray
  - GSE83452 NAFLD spectrum (n = 231), microarray
- Density plot marking the uploaded sample (downloadable as SVG)
- Gene-coverage report (which of the 160 score genes were found)
- Downloadable result CSV

## Scoring details (important for reproducibility)

ssGSEA ranks genes within each sample, but (i) filters genes that are constant
across the input matrix and (ii) divides scores by the global enrichment-score
range computed across all samples. Scoring a sample in isolation is therefore
degenerate. This app scores each upload **jointly with a fixed reference
matrix** (GSE193066, 164 samples, 17,231-gene universe, author-normalized
values log2(x+1), duplicate symbols collapsed by highest expression), which
anchors both steps to the reference cohort.

The ssGSEA implementation in `R/rps_score.R` is a pure-R exact replica of the
GSVA v2 algorithm (validated against GSVA 2.0.0: max |difference| = 0 over all
arm scores of the 164-sample reference), so the app only needs CRAN packages.

Validation (see `data/validation_summary.csv`):

| Check | Result |
|---|---|
| GSE193066 batch recompute from raw GCT (164 samples) | max \|ΔRPS\| = 5.6e-16 |
| Each reference sample re-scored via the app path | max \|ΔRPS\| = 5.6e-16 |
| Example file reproduction | exact |
| External cohort (GSE135251, n = 216) scored via app vs published batch | Spearman ρ = 0.9996, median shift −0.0019 |

Five of the 160 score genes (FYB1, GRK3, JAML, MACROH2A2, TCIM) are absent from
the reference universe and are always reported as missing; this matches the
published GSE193066 analysis (155/160 genes scored).

## Run locally

```r
# R >= 4.3; CRAN packages only: shiny, ggplot2
shiny::runApp("rps_calculator")
```

## Deploy to shinyapps.io

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
rsconnect::deployApp("rps_calculator", appName = "rps-calculator")
```

## Files

- `app.R` — Shiny application
- `R/rps_score.R` — scoring functions (mirror the manuscript pipeline)
- `data/score_genes.csv` — 160-gene score definition (gene, component, regr_logFC, arm)
- `data/reference_expr.rds` — fixed GSE193066 reference expression matrix
- `data/reference_rps.csv` — precomputed RPS distributions of the reference panels
- `data/validation_summary.csv` — validation results
- `example_input.csv` — example upload (GSE193066 HUnafld001)
