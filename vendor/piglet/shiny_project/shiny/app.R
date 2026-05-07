# PIgLET Shiny Application
# Main entry point

# Load global configuration
source("global.R")

# Load modules
source("modules/mod_file_upload.R")
source("modules/mod_clustering.R")
source("modules/mod_visualization.R")

# UI Definition
ui <- navbarPage(
  title = "PIgLET - Allele Similarity Clusters",
  id = "main_nav",
  theme = NULL,

  # Home Tab
  tabPanel(
    "Home",
    icon = icon("home"),
    fluidPage(
      br(),
      fluidRow(
        column(
          width = 8, offset = 2,
          div(
            class = "jumbotron text-center",
            h1("PIgLET"),
            h3("Program for Inferring Immunoglobulin Allele Similarity Clusters"),
            hr(),
            p(
              "PIgLET is a suite of computational tools that improves genotype inference ",
              "and downstream AIRR-seq data analysis. The package provides:"
            ),
            tags$ul(
              class = "list-unstyled",
              tags$li(icon("sitemap"), " Allele Similarity Clusters (ASC) for reducing allele ambiguity"),
              tags$li(icon("dna"), " Support for V, D, and J segments"),
              tags$li(icon("project-diagram"), " Hierarchical and Leiden community detection"),
              tags$li(icon("chart-line"), " Allele-based genotype inference")
            ),
            hr(),
            p("Select 'ASC Inference' above to begin clustering your germline set."),
            br(),
            tags$a(
              href = "https://doi.org/10.1093/nar/gkad603",
              target = "_blank",
              class = "btn btn-info",
              icon("book"), " Read the Paper"
            ),
            tags$a(
              href = "https://github.com/yaarilab/piglet",
              target = "_blank",
              class = "btn btn-secondary",
              icon("github"), " GitHub"
            )
          )
        )
      )
    )
  ),

  # ASC Inference Tab
  tabPanel(
    "ASC Inference",
    icon = icon("sitemap"),
    fluidPage(
      br(),
      fluidRow(
        # Left sidebar - Controls
        column(
          width = 3,
          fileUploadUI("upload"),
          clusteringUI("clustering")
        ),

        # Main content - Results
        column(
          width = 9,
          visualizationUI("visualization")
        )
      )
    )
  ),

  # Documentation Tab
  tabPanel(
    "Documentation",
    icon = icon("book"),
    fluidPage(
      br(),
      fluidRow(
        column(
          width = 10, offset = 1,
          h2("User Guide"),
          hr(),

          h3("Getting Started"),
          p("1. Upload a FASTA file containing your germline sequences, or use the default human V set."),
          p("2. Configure sequence processing options (trimming and masking)."),
          p("3. Choose your clustering method and parameters."),
          p("4. Click 'Run Clustering' to generate allele similarity clusters."),

          hr(),

          h3("Clustering Methods"),

          h4("Hierarchical Clustering (Default)"),
          p("Traditional approach using distance-based hierarchical clustering with two thresholds:"),
          tags$ul(
            tags$li(strong("Family Threshold (75%)"), " - Groups alleles into family-level clusters"),
            tags$li(strong("Allele Cluster Threshold (95%)"), " - Groups alleles into finer clusters")
          ),
          p("Best for V segments with aligned sequences."),

          h4("Leiden Community Detection"),
          p("Modern graph-based community detection algorithm. Benefits:"),
          tags$ul(
            tags$li("Handles variable-length sequences (D and J segments)"),
            tags$li("Automatic resolution optimization using silhouette score"),
            tags$li("Network visualization of cluster relationships")
          ),

          hr(),

          h3("Distance Methods"),
          tags$dl(
            tags$dt("DECIPHER"),
            tags$dd("Default method using DECIPHER package. Best for aligned sequences. Handles gaps properly."),

            tags$dt("Hamming"),
            tags$dd("Simple character-by-character comparison. Sequences padded to equal length."),

            tags$dt("Levenshtein"),
            tags$dd("Edit distance allowing insertions, deletions, substitutions. Good for variable-length sequences.")
          ),

          hr(),

          h3("Output"),
          p("The clustering produces:"),
          tags$ul(
            tags$li(strong("Cluster Table"), " - Allele assignments to clusters"),
            tags$li(strong("Renamed Germline Set"), " - FASTA with ASC nomenclature"),
            tags$li(strong("Visualizations"), " - Circular dendrogram or network plot")
          ),

          hr(),

          h3("Reference"),
          p("Peres A, Lees W, Yan WY, et al. IGHV allele similarity clustering improves genotype inference from AIRR-seq data. Nucleic Acids Research. 2023."),
          tags$a(href = "https://doi.org/10.1093/nar/gkad603", target = "_blank", "DOI: 10.1093/nar/gkad603")
        )
      )
    )
  )
)

# Server Definition
server <- function(input, output, session) {

  # Initialize modules
  upload_data <- fileUploadServer("upload")
  clustering_data <- clusteringServer("clustering", upload_data)
  visualizationServer("visualization", clustering_data, upload_data)

}

# Run the application
shinyApp(ui = ui, server = server)
