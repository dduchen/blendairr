# Global configuration for PIgLET Shiny App
# Loaded once when the app starts

library(shiny)
library(piglet)
library(ggplot2)
library(data.table)
library(DT)

# Load optional packages
have_ggtree <- requireNamespace("ggtree", quietly = TRUE)
have_msaR <- requireNamespace("msaR", quietly = TRUE)
have_ggmsa <- requireNamespace("ggmsa", quietly = TRUE)
have_msa <- requireNamespace("msa", quietly = TRUE)

# Default theme for plots
asc_theme <- function(axis.text.y = 20,
                      strip.text = 40,
                      axis.text.x.bottom = 18,
                      axis.text.x.top = 18,
                      legend.text = 20,
                      legend.title = 20,
                      axis.title.y = 20,
                      axis.title.x = 20) {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      panel.grid = ggplot2::element_line(color = "#cccccc", linewidth = 0.1),
      panel.grid.major = ggplot2::element_line(color = "#cccccc", linewidth = 0.1),
      panel.grid.minor = ggplot2::element_line(color = "#cccccc", linewidth = 0.05),
      axis.ticks = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = axis.text.y, margin = ggplot2::margin(r = 0)),
      strip.text = ggplot2::element_text(size = strip.text, face = "plain", family = ""),
      axis.text.x.bottom = ggplot2::element_text(size = axis.text.x.bottom, angle = 0, vjust = 0.5, hjust = 0.5),
      axis.text.x.top = ggplot2::element_text(size = axis.text.x.top, angle = 0, vjust = 0.5, hjust = 0.5),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(),
      legend.text = ggplot2::element_text(size = legend.text),
      legend.title = ggplot2::element_text(size = legend.title),
      axis.title.y = ggplot2::element_text(size = axis.title.y),
      axis.title.x = ggplot2::element_text(size = axis.title.x)
    )
}

# Locus options for different segments
locus_options <- list(
  "Heavy Chain V" = "IGHV",
  "Heavy Chain D" = "IGHD",
  "Heavy Chain J" = "IGHJ",
  "Kappa Chain V" = "IGKV",
  "Kappa Chain J" = "IGKJ",
  "Lambda Chain V" = "IGLV",
  "Lambda Chain J" = "IGLJ"
)

# Clustering method options
clustering_methods <- list(
  "Hierarchical (default)" = "hierarchical",
  "Leiden Community Detection" = "leiden"
)

# Distance method options
distance_methods <- list(
  "DECIPHER (aligned sequences)" = "decipher",
  "Hamming (equal length)" = "hamming",
  "Levenshtein (variable length)" = "lv"
)

# Helper function
unique_l <- function(l) {
  length(unique(l))
}

# Extract cluster summary from ASC object
extract_clusters <- function(asc) {
  data.frame(
    thresh_fam = asc$threshold[[1]],
    fam_clust = unique_l(asc$alleleClusterTable$family),
    thresh_asc = asc$threshold[[2]],
    asc_clust  = unique_l(asc$alleleClusterTable$allele_cluster)
  )
}
