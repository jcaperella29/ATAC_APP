library(shiny)
library(shinyjs)
library(plotly)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(clusterProfiler)
library(enrichR)
library(DT)

# ==== 🔐 Error Logger ====
log_error <- function(e, context = "unknown") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", timestamp, "] [", context, "] ", conditionMessage(e), "\n")
  write(msg, file = "error_log.txt", append = TRUE)
}

# ==== UI ====
ui <- fluidPage(
  includeCSS("www/fairy_tail.css"),
  useShinyjs(),
  titlePanel("ATAC-seq Peak Annotation & Enrichment Viewer"),
  sidebarLayout(
    sidebarPanel(
      fileInput("peakfile", "Upload MACS2 narrowPeak File", accept = ".narrowPeak"),
      selectInput("enrichr_db", "Select Enrichment Database:",
                  choices = c("GO_Biological_Process_2023",
                              "KEGG_2021_Human",
                              "Reactome_2022"),
                  selected = "GO_Biological_Process_2023"),
      actionButton("run_annotation", "Run Peak Annotation"),
      actionButton("run_enrichment", "Run Enrichment Analysis"),
      div(id = "status_msg", style = "margin-top:10px; font-weight:bold; color:#2c3e50;"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Peak Annotation Table", DTOutput("peak_table"), downloadButton("download_peak_table", "Download CSV")),
        tabPanel("Annotation Pie Chart", plotOutput("pie_plot")),
        tabPanel("Enrichment Table", DTOutput("enrich_table"), downloadButton("download_enrich_table", "Download CSV")),
        tabPanel("Enrichment Bar Plot", plotlyOutput("bar_plot")),
        tabPanel("README", 
                 div(
                   includeMarkdown("HOW_TO.md"),
                   style = "
      max-height: 500px;
      overflow-y: auto;
      background-color: #f9f9f9;
      border-radius: 8px;
      padding: 20px;
      box-shadow: inset 0 0 10px rgba(0,0,0,0.1);
    "
                 )
        )
        
        
      )
    )
  )
)

# ==== Server ====
server <- function(input, output, session) {
  peak_data <- reactiveVal()
  annotation_data <- reactiveVal()
  enrichment_data <- reactiveVal()
  
  showStatus <- function(message) {
    showNotification(message, type = "message", duration = 5)
  }
  
  observeEvent(input$run_annotation, {
    tryCatch({
      req(input$peakfile)
      showStatus("Peak annotation started...")
      peak <- readPeakFile(input$peakfile$datapath)
      peak_data(peak)
      annotation <- annotatePeak(peak,
                                 tssRegion = c(-3000, 3000),
                                 TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                 annoDb = "org.Hs.eg.db")
      annotation_data(annotation)
      showStatus("Peak annotation complete.")
    }, error = function(e) {
      log_error(e, "run_annotation")
      showNotification("❌ Error during peak annotation. Check error_log.txt.", type = "error")
    })
  })
  
  observeEvent(input$run_enrichment, {
    tryCatch({
      req(annotation_data())
      showStatus("Enrichment analysis started...")
      genes <- unique(na.omit(as.data.frame(annotation_data())$SYMBOL))
      req(length(genes) > 0)
      enrichment <- enrichr(genes, input$enrichr_db)
      enrichment_data(enrichment)
      showStatus("Enrichment analysis complete.")
    }, error = function(e) {
      log_error(e, "run_enrichment")
      showNotification("❌ Error during enrichment analysis. Check error_log.txt.", type = "error")
    })
  })
  
  output$peak_table <- renderDT({
    tryCatch({
      req(annotation_data())
      DT::datatable(as.data.frame(annotation_data()))
    }, error = function(e) {
      log_error(e, "render_peak_table")
      NULL
    })
  })
  
  output$download_peak_table <- downloadHandler(
    filename = function() { "peak_annotations.csv" },
    content = function(file) {
      tryCatch({
        write.csv(as.data.frame(annotation_data()), file, row.names = FALSE)
      }, error = function(e) {
        log_error(e, "download_peak_table")
      })
    }
  )
  
  output$pie_plot <- renderPlot({
    tryCatch({
      req(annotation_data())
      plotAnnoPie(annotation_data())
    }, error = function(e) {
      log_error(e, "pie_plot")
      NULL
    })
  })
  
  output$enrich_table <- renderDT({
    tryCatch({
      req(enrichment_data())
      DT::datatable(enrichment_data()[[1]])
    }, error = function(e) {
      log_error(e, "render_enrich_table")
      NULL
    })
  })
  
  output$download_enrich_table <- downloadHandler(
    filename = function() { paste0("enrichment_", input$enrichr_db, ".csv") },
    content = function(file) {
      tryCatch({
        write.csv(enrichment_data()[[1]], file, row.names = FALSE)
      }, error = function(e) {
        log_error(e, "download_enrich_table")
      })
    }
  )
  
  output$bar_plot <- renderPlotly({
    tryCatch({
      req(enrichment_data())
      df <- enrichment_data()[[1]]
      df <- head(df[order(df$Adjusted.P.value), ], 10)
      df <- na.omit(df)
      plot_ly(df, x = ~reorder(Term, -log10(Adjusted.P.value)),
              y = ~-log10(Adjusted.P.value),
              type = 'bar', orientation = 'v') %>%
        layout(title = paste("Top Enriched Terms:", input$enrichr_db),
               xaxis = list(title = "Term", tickangle = -45),
               yaxis = list(title = "-log10(Adjusted P-value)"))
    }, error = function(e) {
      log_error(e, "bar_plot")
      NULL
    })
  })
}

# ==== Run the app ====
shinyApp(ui, server)
