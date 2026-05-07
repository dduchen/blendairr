# Visualization Module for PIgLET Shiny App

#' Visualization Module UI
#'
#' @param id Module namespace ID
visualizationUI <- function(id) {
  ns <- NS(id)

  tabsetPanel(
    id = ns("viz_tabs"),
    type = "tabs",

    tabPanel(
      "Cluster Table",
      icon = icon("table"),
      br(),
      DT::dataTableOutput(ns("cluster_table")),
      br(),
      downloadButton(ns("download_table"), "Download CSV", class = "btn-sm"),
      downloadButton(ns("download_fasta"), "Download FASTA", class = "btn-sm")
    ),

    tabPanel(
      "Circular Dendrogram",
      icon = icon("circle-notch"),
      br(),
      plotOutput(ns("circular_plot"), height = "800px"),
      br(),
      downloadButton(ns("download_circular"), "Download Plot", class = "btn-sm")
    ),

    tabPanel(
      "Network View",
      icon = icon("project-diagram"),
      br(),
      conditionalPanel(
        condition = sprintf("output['%s'] == 'leiden'", ns("has_network")),
        fluidRow(
          column(4,
                 selectInput(ns("network_layout"), "Layout",
                             choices = c("Fruchterman-Reingold" = "fr",
                                         "Kamada-Kawai" = "kk",
                                         "Circle" = "circle"))),
          column(4,
                 selectInput(ns("node_color"), "Color by",
                             choices = c("Cluster" = "cluster",
                                         "Family" = "family"))),
          column(4,
                 checkboxInput(ns("show_labels"), "Show labels", value = TRUE))
        ),
        plotOutput(ns("network_plot"), height = "600px"),
        br(),
        downloadButton(ns("download_network"), "Download Plot", class = "btn-sm")
      ),
      conditionalPanel(
        condition = sprintf("output['%s'] != 'leiden'", ns("has_network")),
        div(class = "alert alert-info",
            icon("info-circle"),
            "Network view is only available for Leiden clustering. ",
            "Select 'Leiden Community Detection' in Clustering Settings.")
      )
    ),

    tabPanel(
      "Summary Statistics",
      icon = icon("chart-bar"),
      br(),
      fluidRow(
        column(6,
               h4("Cluster Summary"),
               verbatimTextOutput(ns("summary_text"))
        ),
        column(6,
               h4("Cluster Size Distribution"),
               plotOutput(ns("size_distribution"), height = "300px")
        )
      ),
      hr(),
      h4("Silhouette Analysis"),
      conditionalPanel(
        condition = sprintf("output['%s'] == 'leiden'", ns("has_network")),
        plotOutput(ns("silhouette_plot"), height = "300px")
      ),
      conditionalPanel(
        condition = sprintf("output['%s'] != 'leiden'", ns("has_network")),
        p("Silhouette analysis is available for Leiden clustering.")
      )
    )
  )
}

#' Visualization Module Server
#'
#' @param id Module namespace ID
#' @param clustering_data Reactive list from clustering module
#' @param germline_data Reactive list from file upload module
visualizationServer <- function(id, clustering_data, germline_data) {
  moduleServer(id, function(input, output, session) {

    # Output to check if network is available
    output$has_network <- reactive({
      req(clustering_data$method())
      clustering_data$method()
    })
    outputOptions(output, "has_network", suspendWhenHidden = FALSE)

    # Cluster table
    output$cluster_table <- DT::renderDataTable({
      req(clustering_data$result())
      DT::datatable(
        clustering_data$result()$alleleClusterTable,
        options = list(
          pageLength = 25,
          scrollX = TRUE
        ),
        filter = "top"
      )
    })

    # Circular dendrogram plot
    output$circular_plot <- renderPlot({
      req(clustering_data$result())
      result <- clustering_data$result()

      if (result$clusteringMethod == "hierarchical") {
        plot(result)
      } else {
        # For Leiden, create a simple cluster visualization
        piglet::plotCommunityNetwork(result,
                                     layout = "circle",
                                     node_color = "cluster",
                                     show_labels = TRUE)
      }
    })

    # Network plot (Leiden only)
    output$network_plot <- renderPlot({
      req(clustering_data$result())
      result <- clustering_data$result()

      if (result$clusteringMethod == "leiden" && !is.null(result$graphObject)) {
        piglet::plotCommunityNetwork(
          result,
          layout = input$network_layout,
          node_color = input$node_color,
          show_labels = input$show_labels
        )
      }
    })

    # Summary text
    output$summary_text <- renderPrint({
      req(clustering_data$result())
      result <- clustering_data$result()
      summ <- summary(result)

      cat("Locus:", summ$locus, "\n")
      cat("Clustering Method:", summ$clustering_method, "\n")
      cat("Number of Alleles:", summ$n_alleles, "\n")
      cat("Number of Clusters:", summ$n_clusters, "\n")
      cat("Number of Families:", summ$n_families, "\n")

      if (summ$clustering_method == "hierarchical") {
        cat("\nThresholds:\n")
        cat("  Family:", summ$family_threshold, "\n")
        cat("  Cluster:", summ$cluster_threshold, "\n")
      } else {
        cat("\nResolution:", summ$resolution_parameter, "\n")
        if (!is.na(summ$silhouette_score)) {
          cat("Silhouette Score:", round(summ$silhouette_score, 4), "\n")
        }
      }
    })

    # Cluster size distribution
    output$size_distribution <- renderPlot({
      req(clustering_data$result())
      result <- clustering_data$result()

      cluster_sizes <- table(result$alleleClusterTable$allele_cluster)
      size_df <- data.frame(
        Cluster = names(cluster_sizes),
        Size = as.integer(cluster_sizes)
      )

      ggplot(size_df, aes(x = Size)) +
        geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
        labs(x = "Cluster Size (# alleles)", y = "Count",
             title = "Distribution of Cluster Sizes") +
        theme_minimal()
    })

    # Download handlers
    output$download_table <- downloadHandler(
      filename = function() {
        paste0("asc_clusters_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(clustering_data$result()$alleleClusterTable, file, row.names = FALSE)
      }
    )

    output$download_fasta <- downloadHandler(
      filename = function() {
        paste0("asc_germline_", Sys.Date(), ".fasta")
      },
      content = function(file) {
        result <- clustering_data$result()
        seqs <- result$alleleClusterSet
        lines <- paste0(">", names(seqs), "\n", seqs)
        writeLines(lines, file)
      }
    )

    output$download_circular <- downloadHandler(
      filename = function() {
        paste0("asc_plot_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        pdf(file, width = 12, height = 12)
        result <- clustering_data$result()
        if (result$clusteringMethod == "hierarchical") {
          plot(result)
        } else {
          print(piglet::plotCommunityNetwork(result, layout = "circle"))
        }
        dev.off()
      }
    )

    output$download_network <- downloadHandler(
      filename = function() {
        paste0("asc_network_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        result <- clustering_data$result()
        p <- piglet::plotCommunityNetwork(
          result,
          layout = input$network_layout,
          node_color = input$node_color,
          show_labels = input$show_labels
        )
        ggsave(file, p, width = 10, height = 8)
      }
    )
  })
}
