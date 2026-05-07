# Contents of app_server.R
asc_theme <- function(axis.text.y = 20,
                      strip.text = 40,
                      axis.text.x.bottom = 18,
                      axis.text.x.top = 18,
                      legend.text = 20,
                      legend.title = 20,
                      axis.title.y = 20,
                      axis.title.x = 20) {
  ggplot2::theme_minimal() +
    theme(
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.grid = element_line(color = "#cccccc", linewidth = 0.1),
      panel.grid.major = element_line(color = "#cccccc", linewidth = 0.1),
      panel.grid.minor = element_line(color = "#cccccc", linewidth = 0.05),
      axis.ticks = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      plot.background = element_blank(),
      axis.text.y = element_text(size = axis.text.y, margin = margin(r = 0)),
      strip.text = element_text(
        size = strip.text,
        face = "plain",
        family = ""
      ),
      axis.text.x.bottom = element_text(
        size = axis.text.x.bottom,
        angle = 0,
        vjust = 0.5,
        hjust = 0.5
      ),
      axis.text.x.top = element_text(
        size = axis.text.x.top,
        angle = 0,
        vjust = 0.5,
        hjust = 0.5
      ),
      strip.placement = "outside",
      strip.background = element_rect(),
      legend.text = element_text(size = legend.text),
      legend.title = element_text(size = legend.title),
      axis.title.y = element_text(size = axis.title.y),
      axis.title.x = element_text(size = axis.title.x)
    )
}

unique_l <- function(l) {
  length(unique(l))
}

extract_clusters <- function(asc) {
  data.frame(
    thresh_fam = asc$threshold[[1]],
    fam_clust = unique_l(asc$alleleClusterTable$family),
    thresh_asc = asc$threshold[[2]],
    asc_clust  = unique_l(asc$alleleClusterTable$allele_cluster)
  )
}

genes_to_cluster <- function(asc) {
  asc$alleleClusterTable$iuis_gene <-
    alakazam::getGene(asc$alleleClusterTable$iuis_allele,
                      collapse = T,
                      strip_d = F)
  asc$alleleClusterTable %>% group_by(family) %>%
    summarise(genes = sum(unique_l(iuis_gene) > 1),
              thresh = unique(asc$threshold[[1]])) %>%
    group_by(thresh) %>%
    summarise(
      clust_gene = sum(genes) / unique_l(family),
      num_clust = unique_l(Family),
      num_shared = sum(genes)
    )
}

asc_thresholds <- function(germline_set) {
  thresholds <- c(0:100)
  
  fam_asc <-
    lapply(thresholds,
           function(i) {
             piglet::inferAlleleClusters(
               germline_set,
               family_threshold = i,
               allele_cluster_threshold = i
             )
           })
  
  fam_asc_tab <- rbindlist(lapply(fam_asc, extract_clusters))
  
  gene_clust_tab <- rbindlist(lapply(fam_asc, genes_to_cluster))
  
  asc_to_iuis <- sapply(fam_asc, function(x) {
    iuis_g <- getGene(x$alleleClusterTable$iuis_allele, strip_d = F)
    asc_g <-
      x$alleleClusterTable %>% mutate(g = getGene(iuis_allele, strip_d = F),
                                      asc = getGene(new_allele, strip_d = F)) %>%
      dplyr::group_by(g) %>%
      dplyr::summarise(asc = unique_l(asc)) %>%
      filter(asc > 1) %>% pull(asc) %>% length()
    asc_g / unique_l(iuis_g)
  })
  
  tab <-
    data.frame(
      thresh = fam_asc_tab$thresh_fam,
      clust_alleles = fam_asc_tab$fam_clust / length(germline_set),
      gene_clust = gene_clust_tab$clust_gene,
      iuis_asc = asc_to_iuis
    )
  
  return(tab)
}

server <- function(input, output, session) {
  germline_set <- reactiveVal(NULL)
  inferred_asc <- reactiveVal(NULL)
  
  observeEvent(input$file, {
    raw_sequences <- tigger::readIgFasta(input$file$datapath)
    
    germline_set(unlist(raw_sequences))
    ## check if alined. if not align.
    check <- all(grepl("[.]", raw_sequences))
    if (!check) {
      print("aligning sequences")
      germline_set(piglet:::alignSeqs(raw_sequences))
    }
    
    
    
    max_val <- max(nchar(germline_set()))
    min_val <- min(nchar(germline_set())) - 50
    
    # Set the range of trim_3prime_side slider
    updateSliderInput(
      session,
      "trim_3prime_side",
      min = min_val,
      max = max_val,
      value = max_val - 2
    )
  })
  
  observeEvent(input$use_default_germline, {
    germline_set(tigger::readIgFasta(system.file("data", "HVGERM.fasta", package =
                                                   "piglet")))
    
    # Set the range of trim_3prime_side slider
    max_val <- max(nchar(germline_set()))
    min_val <- min(nchar(germline_set())) - 50
    
    # Set the range of trim_3prime_side slider
    updateSliderInput(
      session,
      "trim_3prime_side",
      min = min_val,
      max = max_val,
      value = max_val - 2
    )
  })
  
  # observe({
  #   # Enforce condition for trim_3prime_side and mask_5prime_side visibility
  #   if (!is.null(input$file) ||
  #       input$use_default_germline > 0) {
  #     shinyjs::enable("trim_3prime_side")
  #     shinyjs::enable("mask_5prime_side")
  #     shinyjs::enable("family_threshold")
  #     shinyjs::enable("allele_cluster_threshold")
  #   } else {
  #     shinyjs::disable("trim_3prime_side")
  #     shinyjs::disable("mask_5prime_side")
  #     shinyjs::disable("family_threshold")
  #     shinyjs::disable("allele_cluster_threshold")
  #   }
  # })
  
  observe({
    # Enforce condition for trim_3prime_side
    if (!is.null(germline_set())) {
      if (input$mask_5prime_side > input$trim_3prime_side) {
        updateSliderInput(session, "trim_3prime_side", value = input$mask_5prime_side)
      }
      
      # Enforce condition for mask_5prime_side
      if (input$trim_3prime_side < input$mask_5prime_side) {
        updateSliderInput(session, "mask_5prime_side", value = input$trim_3prime_side)
      }
    }
  })
  
  ## MSA plot ----------------------------------------------
  observe({
    if (!is.null(germline_set())) {
      raw_sequences <- isolate(germline_set())
      
      ## If seqs are not the same length. Add alignment dots.
      max_length <- max(nchar(raw_sequences))
      padded_sequences <-
        unlist(sapply(raw_sequences, function(seq) {
          if (nchar(seq) < max_length) {
            padded_seq <-
              paste0(seq, paste0(rep(".", max_length - nchar(seq)), collapse = ""), collapse = "")
            padded_seq
          } else {
            seq
          }
        }))
      padded_sequences <-
        Biostrings::AAStringSet(gsub("[.]", "-", padded_sequences))
      
      output$msa <- renderMsaR({
        msaR(padded_sequences,
             menu = FALSE,
             seqlogo = TRUE)
      })
      
      output$download_msa_pdf <- downloadHandler(
        filename = function() {
          paste("msa_germline_set_", Sys.Date(), ".pdf", sep = "")
        },
        content = function(file) {
          msa::msaPrettyPrint(msa::msa(padded_sequences) 
                              , file = 'myreport.pdf'
                              , output="pdf"
                              , showNames="left"
                              , showLogo="top"
                              , consensusColor="BlueRed"
                              , logoColors="accessible area"
                              , askForOverwrite=FALSE)
          file.rename('myreport.pdf', file)
        },
        contentType = 'application/pdf'
      )
      
      output$download_msa_png <- downloadHandler(
        filename = function() {
          paste("msa_germline_set_", Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
          p <- ggmsa(Biostrings::DNAStringSet(germline_set()), 
                char_width = 0.5, seq_name = T, custom_color = data.frame(
            names = c("A", "T", "C", "G", "N", "-"),
            color = c(
              "#ff6d6d",
              "#769dcc",
              "#f2be3c",
              "#74ce98",
              "#b8b8b8",
              "#ffffff"
            )
          )) + geom_seqlogo() + geom_msaBar()
          ggsave(file, plot = p, device = "png")
        },
        contentType = 'application/pdf'
      )
    }
  })
  
  ## 3' trimming selection ---------
  observe({
    if (!is.null(germline_set())) {
      ## An assist to the 3' trimming selection
      
      
      # upper bound
      germline_set_ends <- nchar(isolate(germline_set()))
      max_val <- max(germline_set_ends)
      selected_val <- input$trim_3prime_side
      # get ecdf
      Fn_end <- ecdf(germline_set_ends)
      output$trim3 <- renderPlot({
        ggplot() +
          geom_point(
            mapping = aes(x, y),
            data = data.frame(
              x = c(selected_val, max_val),
              y = c(Fn_end(selected_val), Fn_end(max_val))
            ),
            size = 3
          ) +
          stat_ecdf(
            mapping = aes(x),
            data = data.frame(x = germline_set_ends),
            geom = "step"
          ) +
          geom_vline(
            xintercept = c(selected_val, max_val),
            color = "red",
            linetype = "dashed"
          ) +
          geom_hline(
            yintercept = c(Fn_end(selected_val), Fn_end(max_val)),
            color = "gray70",
            linetype = "dashed"
          ) +
          asc_theme() +
          labs(x = "Allele sequence most 3' position",
               y = "ecdf(position)")
      })
      
    }
  })
  
  
  ## ASC thresold inference ---------
  observe({
    if (!is.null(germline_set())) {
      ## the asc will be inferred to all similarity thresholds,
      ## and a range of 3' trimming. From max to 2 nucliotides below minimum.
      
      tab <- asc_thresholds(germline_set())
      
      output$download_asc_threshols_table <- downloadHandler(
        filename = function() {
          paste("asc_thresholds_table_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(tab, file)
        }
      )
      
      #f1 <- approxfun(tab$gene_clust - tab$clust_alleles, tab$thresh, rule=2)
      rec_thresh_val <-
        tab$thresh[which.min(tab$gene_clust + tab$iuis_asc)] #round(f1(0))
      intercept_vlines <- c(rec_thresh_val,
                            input$allele_cluster_threshold,
                            input$family_threshold)
      
      data_vline <- data.frame(
        line =  c(
          paste0("recommended\nasc\nthreshold (", rec_thresh_val, ")"),
          "selected\nasc\nthreshold",
          "selected\nfamily\nthreshold"
        ),
        value  = intercept_vlines,
        color = c("firebrick",
                  "darkgreen",
                  "deeppink4")
      )
      
      output$asc_threshold <- renderPlot({
        ggplot(tab, aes(thresh)) +
          geom_rect(
            aes(
              xmin = rec_thresh_val - 2,
              xmax = rec_thresh_val + 2,
              ymin = -Inf,
              ymax = Inf
            ) ,
            fill = "gray90"
          ) +
          geom_line(aes(y = iuis_asc), color = "#018571") +
          geom_point(aes(y = iuis_asc)) +
          geom_line(aes(y = gene_clust), color = "#a6611a") +
          geom_point(aes(y = gene_clust)) +
          geom_line(aes(y = gene_clust + iuis_asc), color = "black") +
          geom_point(aes(y = gene_clust + iuis_asc)) +
          scale_y_continuous(
            name = "Num. clusters sharing IUIS \n genes over num. clusters",
            sec.axis = sec_axis(~ . * 1, name = "Num. genes spanning more than \n one cluster over the num. genes")
          ) +
          geom_vline(
            data = data_vline,
            mapping = aes(xintercept = value, colour = line),
            linetype = "dashed"
          ) +
          scale_color_manual(name = "",
                             values = setNames(alpha(data_vline$color, alpha = 1), data_vline$line)) +
          geom_hline(
            yintercept = tab$iuis_asc[tab$thresh  %in% intercept_vlines],
            color = "#5ab4ac",
            linetype = "dashed"
          ) +
          geom_hline(
            yintercept = tab$gene_clust[tab$thresh  %in% intercept_vlines],
            color = "#dfc27d",
            linetype = "dashed"
          ) +
          labs(x = "Similarity threshold") +
          asc_theme(
            axis.text.y = 16,
            axis.title.y = 16,
            axis.text.x.bottom = 16,
            axis.title.x = 16,
            legend.text = 14
          ) +
          theme(
            axis.text.y.right = element_text(color = "#018571"),
            axis.text.y.left = element_text(color = "#a6611a"),
            legend.position = "bottom"
          )
      })
      
    }
  })
  
  ## ASC inference ------------
  ## ASC inference ---------
  observe({
    if (!is.null(germline_set())) {
      inferred_asc_result <- piglet::inferAlleleClusters(
        germline_set(),
        family_threshold = input$family_threshold,
        allele_cluster_threshold = input$allele_cluster_threshold,
        trim_3prime_side = input$trim_3prime_side,
        mask_5prime_side = input$mask_5prime_side
      )
      inferred_asc(inferred_asc_result)
      output$asc_table <- DT::renderDataTable(inferred_asc()$alleleClusterTable)
      
      output$asc_plot <- renderPlot({
        plot(inferred_asc())
      },)
      
    }
  })
  
  ## Download buttons ------------
  output$download_asc_table <- downloadHandler(
    filename = function() {
      paste("asc_table_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(inferred_asc()$alleleClusterTable, file)
    }
  )
  
  output$download_fasta <- downloadHandler(
    filename = function() {
      paste("germline_set_", Sys.Date(), ".fasta", sep = "")
    },
    content = function(file) {
      writeLines(germlineASC(inferred_asc()$alleleClusterTable, germline_set()),
                 file)
    }
  )
  
  
}
