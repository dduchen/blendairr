# File Upload Module for PIgLET Shiny App

#' File Upload Module UI
#'
#' @param id Module namespace ID
fileUploadUI <- function(id) {
  ns <- NS(id)

  wellPanel(
    h4("Input Data"),
    fluidRow(
      column(width = 6,
             fileInput(ns("file"), "Upload FASTA file",
                       accept = c(".fasta", ".fa", ".fna"))),
      column(width = 6,
             actionButton(ns("use_default"), "Use Default Human V Set",
                          class = "btn-primary btn-sm",
                          style = "margin-top: 25px;"))
    ),
    hr(),
    h4("Sequence Processing"),
    sliderInput(
      ns("trim_3prime"),
      "Trim 3' Side (position)",
      min = 0,
      max = 500,
      value = 318
    ),
    sliderInput(
      ns("mask_5prime"),
      "Mask 5' Side (nucleotides)",
      min = 0,
      max = 200,
      value = 0
    )
  )
}

#' File Upload Module Server
#'
#' @param id Module namespace ID
#' @return List of reactive values: germline_set, trim_3prime, mask_5prime
fileUploadServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Reactive value to store germline set
    germline_set <- reactiveVal(NULL)

    # Handle file upload
    observeEvent(input$file, {
      req(input$file)

      tryCatch({
        raw_sequences <- tigger::readIgFasta(input$file$datapath)
        germline_set(unlist(raw_sequences))

        # Check if aligned, if not show warning
        check <- all(grepl("[.]", raw_sequences))
        if (!check) {
          showNotification(
            "Sequences appear unaligned. For best results, use IMGT-gapped sequences.",
            type = "warning",
            duration = 5
          )
        }

        # Update slider range based on sequence lengths
        max_val <- max(nchar(germline_set()))
        min_val <- max(0, min(nchar(germline_set())) - 50)

        updateSliderInput(session, "trim_3prime",
                          min = min_val, max = max_val,
                          value = min(318, max_val))

        showNotification(
          paste("Loaded", length(raw_sequences), "sequences"),
          type = "message"
        )

      }, error = function(e) {
        showNotification(
          paste("Error loading file:", e$message),
          type = "error"
        )
      })
    })

    # Handle default germline button
    observeEvent(input$use_default, {
      tryCatch({
        default_file <- system.file("data", "HVGERM.fasta", package = "piglet")
        if (default_file == "") {
          # Try loading from package data
          data("HVGERM", package = "piglet", envir = environment())
          germline_set(HVGERM)
        } else {
          germline_set(tigger::readIgFasta(default_file))
        }

        max_val <- max(nchar(germline_set()))
        min_val <- max(0, min(nchar(germline_set())) - 50)

        updateSliderInput(session, "trim_3prime",
                          min = min_val, max = max_val,
                          value = min(318, max_val))

        showNotification(
          paste("Loaded default HVGERM set with", length(germline_set()), "sequences"),
          type = "message"
        )

      }, error = function(e) {
        showNotification(
          paste("Error loading default data:", e$message),
          type = "error"
        )
      })
    })

    # Enforce constraint: mask <= trim
    observe({
      req(germline_set())
      if (input$mask_5prime > input$trim_3prime) {
        updateSliderInput(session, "mask_5prime", value = input$trim_3prime)
      }
    })

    # Return reactive values
    list(
      germline_set = germline_set,
      trim_3prime = reactive(input$trim_3prime),
      mask_5prime = reactive(input$mask_5prime)
    )
  })
}
