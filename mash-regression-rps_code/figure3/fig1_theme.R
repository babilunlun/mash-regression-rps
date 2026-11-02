# Figure 1 shared theme — per ipf-lf-paper-rules skill (visual conventions only)
suppressMessages(library(ggplot2))

palette8 <- c(
  "#DCA7EB", "#B1C4D7", "#BFDA8E", "#E2EA4D",
  "#BEF3C5", "#C8CACB", "#C7B0C7", "#A7C0DF"
)
ALPHA <- 0.75

theme_ipf_lf <- function(base = 12) {
  theme_classic(base_size = base) +
    theme(
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      axis.text          = element_text(color = "black", size = base - 1),
      axis.title         = element_text(color = "black", size = base, face = "bold"),
      axis.line          = element_line(color = "black", linewidth = 0.5),
      axis.ticks         = element_line(color = "black", linewidth = 0.4),
      plot.title         = element_text(face = "bold", size = base + 1, hjust = 0.5, color = "black"),
      plot.subtitle      = element_text(size = base - 3, hjust = 0.5, color = "#666666"),
      legend.background  = element_rect(fill = "white", color = NA),
      panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      strip.text         = element_text(face = "bold", size = base, color = "black"),
      strip.background   = element_rect(fill = "#F5F5F5", color = "#CCCCCC", linewidth = 0.5),
      plot.margin        = margin(12, 40, 12, 12),
      text               = element_text(family = "Liberation Sans")
    )
}

# Save SVG + PNG pair. width/height default to single-panel 5x4 in.
save_panel <- function(plot, path_stem, width = 5, height = 4) {
  svglite::svglite(paste0(path_stem, ".svg"), width = width, height = height)
  print(plot); dev.off()
  png(paste0(path_stem, ".png"), width = width, height = height,
      units = "in", res = 300, type = "cairo")
  print(plot); dev.off()
  invisible(NULL)
}
