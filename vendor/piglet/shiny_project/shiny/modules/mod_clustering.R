# Clustering Module for PIgLET Shiny App

#' Clustering Module UI
#'
#' @param id Module namespace ID
clusteringUI <- function(id) {
  ns <- NS(id)

  wellPanel(
    h4("Clustering Settings"),

    selectInput(
      ns("locus"),
      "Locus Type",
      choices = locus_options,
      selected = "IGHV"
    ),

    radioButtons(
      ns("method"),
      "Clustering Method",
      choices = clustering_methods,
      selected = "hierarchical"
    ),

    selectInput(
      ns("distance_method"),
      "Distance Method",
      choices = distance_methods,
      selected = "decipher"
    ),

    hr(),

    # Hierarchical clustering parameters
    conditionalPanel(
      condition = sprintf("input['%s'] == 'hierarchical'", ns("method")),
      h5("Hierarchical Clustering"),
      sliderInput(
        ns("family_threshold"),
        "Family Threshold (%)",
        min = 0, max = 100, value = 75
      ),
      sliderInput(
        ns("allele_cluster_threshold"),
        "Allele Cluster Threshold (%)",
        min = 0, max = 100, value = 95
      ),
      selectInput(
        ns("linkage_method"),
        "Linkage Method",
        choices = c("complete", "average", "single", "ward.D2"),
        selected = "complete"
      )
    ),

    # Leiden clustering parameters
    conditionalPanel(
      condition = sprintf("input['%s'] == 'leiden'", ns("method")),
      h5("Leiden Community Detection"),
      numericInput(
        ns("resolution"),
        "Resolution (leave empty for auto)",
        value = NA,
        min = 0.01, max = 10, step = 0.1
      ),
      numericInput(
        ns("target_clusters"),
        "Target # Clusters (for auto resolution)",
        value = NA,
        min = 2, max = 500
      ),
      checkboxInput(
        ns("optimize_silhouette"),
        "Optimize using silhouette score",
        value = TRUE
      ),
      numericInput(
        ns("ncores"),
        "Number of CPU cores",
        value = 1,
        min = 1, max = 16
      )
    ),

    hr(),

    actionButton(ns("run"), "Run Clustering",
                 class = "btn-success btn-block",
                 icon = icon("play"))
  )
}

#' Clustering Module Server
#'
#' @param id Module namespace ID
#' @param germline_data Reactive list from file upload module
#' @return Reactive GermlineCluster object
clusteringServer <- function(id, germline_data) {
  moduleServer(id, function(input, output, session) {

    # Reactive value for results
    asc_result <- reactiveVal(NULL)

    # Run clustering when button clicked
    observeEvent(input$run, {
      req(germline_data$germline_set())

      # Show progress
      withProgress(message = "Running clustering...", value = 0, {

        incProgress(0.1, detail = "Preparing sequences")

        germline_set <- germline_data$germline_set()
        trim_3prime <- germline_data$trim_3prime()
        mask_5prime <- germline_data$mask_5prime()

        tryCatch({
          incProgress(0.3, detail = "Calculating distances")

          # Build parameters
          params <- list(
            germline_set = germline_set,
            locus = input$locus,
            clustering_method = input$method,
            distance_method = input$distance_method,
            trim_3prime_side = if (trim_3prime > 0) trim_3prime else NULL,
            mask_5prime_side = mask_5prime,
            quiet = FALSE
          )

          # Add method-specific parameters
          if (input$method == "hierarchical") {
            params$family_threshold <- input$family_threshold
            params$allele_cluster_threshold <- input$allele_cluster_threshold
            params$cluster_method <- input$linkage_method
          } else {
            params$resolution <- if (is.na(input$resolution)) NULL else input$resolution
            params$target_clusters <- if (is.na(input$target_clusters)) NULL else input$target_clusters
            params$optimize_silhouette <- input$optimize_silhouette
            params$ncores <- input$ncores
          }

          incProgress(0.5, detail = "Running clustering algorithm")

          # Run inference
          result <- do.call(piglet::inferAlleleClusters, params)

          incProgress(0.9, detail = "Finalizing")

          asc_result(result)

          showNotification(
            sprintf("Clustering complete! Found %d clusters with %d alleles",
                    length(unique(result$alleleClusterTable$allele_cluster)),
                    nrow(result$alleleClusterTable)),
            type = "message",
            duration = 5
          )

        }, error = function(e) {
          showNotification(
            paste("Clustering error:", e$message),
            type = "error",
            duration = 10
          )
        })
      })
    })

    # Return reactive result
    list(
      result = asc_result,
      method = reactive(input$method),
      locus = reactive(input$locus)
    )
  })
}
