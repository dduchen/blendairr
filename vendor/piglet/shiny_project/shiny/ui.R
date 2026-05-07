# Contents of app_ui.R

ui <- fluidPage(navbarPage(
  "ASC Inference App",
  tabPanel("Home"),
  tabPanel("ASC Inference",
           fluidRow(
             column(width = 3,
                    wellPanel(
                      fluidRow(column(width = 4,
                                      h5("Choose a file")),
                               column(
                                 width = 6,
                                 actionButton("use_default_germline", "Default human V set", class = "btn-sm")
                               )),
                      fileInput("file", ""),
                      #conditionalPanel(
                        #condition = "input.file || input.use_default_germline",
                        sliderInput(
                          "trim_3prime_side",
                          "Trim 3' Side",
                          min = 0,
                          max = 500,
                          value = 318
                        ),
                        sliderInput(
                          "mask_5prime_side",
                          "Mask 5' Side",
                          min = 0,
                          max = 200,
                          value = 0
                        ),
                        sliderInput(
                          "family_threshold",
                          "Family Threshold",
                          min = 0,
                          max = 100,
                          value = 75
                        ),
                        sliderInput(
                          "allele_cluster_threshold",
                          "Allele Cluster Threshold",
                          min = 0,
                          max = 100,
                          value = 95
                        )
                      #)
                    )),
             column(
               width = 9,
               tabsetPanel(
                 tabPanel("Germline Set",
                          msaROutput("msa"),
                          fluidRow(column(width = 4,
                                          downloadButton("download_msa_pdf", "Download MSA pdf")),
                                   column(
                                     width = 4,
                                     downloadButton("download_msa_png", "Download MSA png")
                                   ))),
                 tabPanel("3' Trimming selection",
                          plotOutput("trim3")),
                 tabPanel("ASC threshold selection",
                          plotOutput("asc_threshold"),
                          downloadButton("download_asc_threshols_table", "Download ASC Table")),
                 tabPanel("Inferred ASC",
                          DT::dataTableOutput("asc_table"),
                          downloadButton("download_asc_table", "Download ASC Table"),
                          downloadButton("download_fasta", "Download Fasta")
                 ),
                 tabPanel("Plot ASC",
                          plotOutput("asc_plot", width = "auto", height = "200vh")
                 )
               )
             )
           ))
))
