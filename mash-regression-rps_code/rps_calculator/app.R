# app.R — Regression Potential Score (RPS) online calculator
# Upload a single-sample gene-expression vector (CSV: gene, expression) and
# obtain the RPS together with its percentile against reference liver cohorts.
#
# RPS = ssGSEA(up arm, 22 genes) - ssGSEA(down arm, 138 genes)
# Arms are the strict cross-species conserved regression core (160 genes).
#
# Uploaded samples are scored JOINTLY with the fixed GSE193066 reference
# matrix (164 samples, 17,231-gene universe). ssGSEA ranks genes within each
# sample but normalizes by the global score range across the input matrix and
# filters genes constant across it; scoring a sample alone is therefore
# degenerate. Anchoring to the fixed reference makes uploaded scores exactly
# comparable to the published cohort values (validated: max |difference| = 0
# for all 164 reference samples). The ssGSEA implementation (R/rps_score.R)
# is a pure-R exact replica of GSVA v2, so the app has no Bioconductor
# dependency.

library(shiny)
library(ggplot2)

source("R/rps_score.R")

arms <- rps_load_genes("data/score_genes.csv")
ref <- read.csv("data/reference_rps.csv", stringsAsFactors = FALSE)
ref_mat <- readRDS("data/reference_expr.rds")
panel_choices <- unique(ref$panel)

theme_set(theme_bw(base_family = "Liberation Sans") +
            theme(text = element_text(family = "Liberation Sans")))

ui <- fluidPage(
  titlePanel("Regression Potential Score (RPS) calculator"),
  p("Upload a normalized liver expression profile (one sample) to compute the RPS, ",
    "a transcriptomic score of the MASH fibrosis regression program. ",
    "Higher RPS indicates a liver expression state closer to fibrosis regression."),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Expression vector (CSV)",
                accept = c(".csv", ".tsv", ".txt")),
      helpText("Two columns: gene symbol and normalized expression ",
               "(e.g. logCPM or log2-normalized intensity). ssGSEA uses ",
               "within-sample ranks, so any monotone transformation of the ",
               "values gives the identical score. An example file is ",
               "bundled with the app."),
      actionButton("example", "Load example input"),
      hr(),
      selectInput("panel", "Reference panel for percentile", choices = panel_choices),
      hr(),
      downloadButton("dl_csv", "Download result (CSV)"),
      downloadButton("dl_svg", "Download plot (SVG)"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Result",
                 br(),
                 fluidRow(
                   column(4, wellPanel(h4("RPS"), h2(textOutput("rps_value")))),
                   column(4, wellPanel(h4("Percentile in reference"),
                                       h2(textOutput("pct_value")))),
                   column(4, wellPanel(h4("Gene coverage"),
                                       h2(textOutput("cov_value"))))
                 ),
                 plotOutput("dens", height = "380px"),
                 br(),
                 verbatimTextOutput("coverage_detail")),
        tabPanel("About",
                 br(),
                 p("The Regression Potential Score contrasts two arms of the ",
                   "cross-species conserved MASH fibrosis regression core: ",
                   "22 genes that increase during regression (up arm) minus ",
                   "138 genes that decrease during regression (down arm), ",
                   "each summarized by single-sample GSEA (ssGSEA)."),
                 p("Because ssGSEA only uses within-sample ranks, the score is ",
                   "portable across platforms (RNA-seq and microarray) and ",
                   "invariant to monotone transformations of expression values."),
                 p("Scoring details: the uploaded vector is de-duplicated by ",
                   "gene symbol (highest value kept), intersected with the ",
                   "17,231-gene reference universe, and scored jointly with ",
                   "the fixed GSE193066 reference matrix (164 samples). This ",
                   "anchors GSVA's constant-gene filter and global score-range ",
                   "normalization to the reference cohort, so uploaded scores ",
                   "are directly comparable to the published values; five of ",
                   "the 160 score genes are absent from the reference universe ",
                   "and are always reported as missing."),
                 p("Reference panels: baseline and follow-up biopsies of a ",
                   "bariatric-surgery cohort with paired histology (GSE193066), ",
                   "and three NAFLD/MASH spectrum cohorts profiled by RNA-seq ",
                   "(GSE135251) or microarray (GSE89632, GSE83452)."),
                 p("Interpretation: a percentile above 50 means the uploaded ",
                   "sample scores higher than half of the chosen reference ",
                   "cohort. Within-patient increases in RPS track fibrosis ",
                   "regression; absolute values are cohort-dependent."))
      ),
      width = 9
    )
  )
)

server <- function(input, output, session) {

  input_data <- reactiveVal(NULL)

  observeEvent(input$file, {
    df <- tryCatch(read.csv(input$file$datapath, stringsAsFactors = FALSE,
                            check.names = FALSE),
                   error = function(e) NULL)
    input_data(df)
  })

  observeEvent(input$example, {
    input_data(read.csv("example_input.csv", stringsAsFactors = FALSE,
                        check.names = FALSE))
  })

  result <- reactive({
    df <- input_data()
    req(df)
    parsed <- rps_parse_input(df)
    validate(need(length(parsed$values) > 100,
                  "Could not parse the input: expected a gene column and a numeric expression column."))
    rps_score_upload(parsed$values, parsed$genes, arms, ref_mat)
  })

  ref_rps <- reactive({
    ref$RPS[ref$panel == input$panel]
  })

  percentile <- reactive({
    r <- result()
    round(100 * ecdf(ref_rps())(r$score$RPS), 1)
  })

  output$rps_value <- renderText({
    sprintf("%.3f", result()$score$RPS)
  })

  output$pct_value <- renderText({
    sprintf("%.1f%%", percentile())
  })

  output$cov_value <- renderText({
    r <- result()
    sprintf("%d / %d", r$n_up + r$n_down,
            length(arms$up) + length(arms$down))
  })

  output$coverage_detail <- renderText({
    r <- result()
    paste0(
      "Up arm found: ", r$n_up, " / ", length(arms$up), "\n",
      "Down arm found: ", r$n_down, " / ", length(arms$down), "\n",
      "Input genes after de-duplication: ", r$n_input_genes, "\n",
      "Genes overlapping the reference universe: ", r$n_universe, " / ",
      nrow(ref_mat), "\n",
      if (length(r$missing_up) > 0)
        paste0("Missing up-arm genes: ", paste(r$missing_up, collapse = ", "), "\n") else "",
      if (length(r$missing_down) > 0)
        paste0("Missing down-arm genes: ", paste(r$missing_down, collapse = ", ")) else ""
    )
  })

  dens_plot <- reactive({
    r <- result()
    rr <- ref_rps()
    x <- r$score$RPS
    ggplot(data.frame(RPS = rr), aes(x = RPS)) +
      geom_density(fill = "grey70", alpha = 0.6, color = "grey30") +
      geom_vline(xintercept = x, color = "#0279EE", linewidth = 1.1,
                 linetype = "dashed") +
      annotate("text", x = x, y = Inf, vjust = 1.5, hjust = -0.1,
               label = sprintf("uploaded sample\nRPS = %.3f\npercentile = %.1f%%",
                               x, percentile()),
               color = "#0279EE", size = 4) +
      labs(x = "RPS", y = "Density",
           title = input$panel) +
      theme(plot.title = element_text(size = 11))
  })

  output$dens <- renderPlot(dens_plot())

  output$dl_csv <- downloadHandler(
    filename = function() "rps_result.csv",
    content = function(file) {
      r <- result()
      out <- data.frame(sample = "uploaded_sample",
                        up = r$score$up, down = r$score$down,
                        RPS = r$score$RPS,
                        reference_panel = input$panel,
                        percentile = percentile(),
                        up_arm_coverage = paste0(r$n_up, "/", length(arms$up)),
                        down_arm_coverage = paste0(r$n_down, "/", length(arms$down)))
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$dl_svg <- downloadHandler(
    filename = function() "rps_density.svg",
    content = function(file) {
      ggsave(file, dens_plot(), device = "svg", width = 7, height = 4.5)
    }
  )
}

shinyApp(ui, server)
