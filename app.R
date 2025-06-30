# app.R
library(shiny); library(shinyjs); library(plotly); library(DT);
 library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db)
library(GenomicRanges);  
library(uwot); library(randomForest); library(pROC)
library(TFBSTools)
library(umap)
library(JASPAR2020)
library(motifmatchr)

library(BSgenome.Hsapiens.UCSC.hg38)
options(shiny.maxRequestSize = 200*1024^2)

# ==== UI ====
ui <- fluidPage(
  includeCSS("www/fairy_tail.css"),
  useShinyjs(),
  titlePanel("JCAP_ATAC_SEQ APP"),
  sidebarLayout(
    sidebarPanel(
      fileInput("zipfile", "Upload ZIP of BED files for consensus peaks", accept = ".zip"),
      actionButton("make_peaks", "Make Consensus Peaks"),
      actionButton("run_annotation", "Run Peak Annotation"),
      actionButton("run_motif", "Run Motif Enrichment"),
      selectInput("jaspar_family", "TF Family (JASPAR)", choices = c("ALL","MAF","TP53","C2H2")),
      fileInput("count_matrix_file", "Upload Count Matrix CSV", accept = ".csv"),
      fileInput("metadata_file", "Upload Metadata CSV", accept = ".csv"),
      actionButton("run_daa_real", "Run DAA on Uploaded Counts"),
      actionButton("run_pca", "Run PCA"),
      actionButton("run_umap", "Run UMAP"),
      actionButton("run_heatmap", "Plot Heatmap"),
      actionButton("run_rf", "Run Random Forest Classifier"),
     
      
      h4("Power Analysis"),
      numericInput("power_effect_size", "Effect Size (log2FC):", value=1, min=0.2, max=3),
      numericInput("power_alpha", "FDR Threshold (α):", value=0.05, min=0.001, max=0.1),
      numericInput("power_replicates", "Replicates per Group:", value=3, min=2, max=20),
      actionButton("run_power", "Run Power Analysis"),
      div(id="status_msg", style="margin-top:10px;font-weight:bold;color:#2c3e50;")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("README", includeMarkdown("HOW_TO.md")),
        tabPanel("Peak Annotation Table", DTOutput("peak_table"), downloadButton("download_peak_table","Download")),
        tabPanel("Annotation Pie Chart", plotOutput("pie_plot")),
        tabPanel("Motif Enrichment Plot", plotlyOutput("motif_enrich_plot")),
        tabPanel("Motif Enrichment Table", DTOutput("motif_enrich_table"), downloadButton("download_motif","Download")),
        tabPanel("Uploaded DAA Results", DTOutput("uploaded_daa_table"), downloadButton("download_uploaded_daa","Download")),
        tabPanel("PCA Plot", plotlyOutput("pca_plot")),
        tabPanel("UMAP Plot", plotlyOutput("umap_plot")),
        tabPanel("Heatmap", plotlyOutput("heatmap_plot")),
        
        tabPanel("RF Metrics", DTOutput("rf_metrics")),
        tabPanel("Feature Importance", plotlyOutput("rf_varimp_plot")),
        tabPanel("ROC Curve", plotlyOutput("roc_plot")),
        tabPanel("Power Plot", plotlyOutput("power_plot")),
        tabPanel("Power Table", DTOutput("power_table"), downloadButton("download_power_table","Download Power Table"))
      )
    )
  )
)

# ==== SERVER ====
server <- function(input, output, session){
  consensus_peaks_rv <- reactiveVal()
  annotation_data <- reactiveVal()
  enrichment_data <- reactiveVal()
  uploaded_counts_rv <- reactiveVal()
  uploaded_daa_rv <- reactiveVal()
  pca_data_rv <- reactiveVal(); umap_data_rv <- reactiveVal()
  motif_enrich_rv <- reactiveVal()
  rf_results_rv <- reactiveVal(); rf_importance_rv <- reactiveVal()
  power_df_rv <- reactiveVal()
  
  showStatus <- function(msg) showNotification(msg, type="message", duration=5)
  log_error <- function(e, ctx) write(paste(Sys.time(), ctx, conditionMessage(e)), file="error_log.txt", append=TRUE)
  
  observeEvent(input$make_peaks, {
    tryCatch({
      req(input$zipfile)
      dir <- tempfile(); unzip(input$zipfile$datapath, exdir = dir)
      beds <- list.files(dir, "\\.(bed|narrowPeak)$", full.names = TRUE, recursive = TRUE)
      
      read_bed_as_granges <- function(file_path) {
        tryCatch({
          df <- read.table(file_path, sep = "\t", header = FALSE)
          if (ncol(df) < 3) stop("File does not have at least 3 columns")
          GRanges(seqnames = df[[1]],
                  ranges = IRanges(start = df[[2]] + 1, end = df[[3]]))
        }, error = function(e) {
          warning(paste("❌ Failed to read BED:", file_path, "Reason:", conditionMessage(e)))
          NULL
        })
      }
      
      seqs <- lapply(beds, read_bed_as_granges)
      seqs <- Filter(Negate(is.null), seqs)
      req(length(seqs) > 0)
      
      consensus_peaks_rv(reduce(do.call(c, seqs)))
      showStatus(paste("✅ Consensus peaks:", length(consensus_peaks_rv())))
      
    }, error = function(e) {
      log_error(e, "make_peaks")
      showNotification("❌ Consensus peak generation failed", type = "error")
    })
  })
  
  # Required libraries
  library(GenomicRanges)
  library(GenomicFeatures)
  library(org.Hs.eg.db)
  library(DT)
  
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes_gr <- genes(txdb)
  tss_gr <- resize(genes_gr, 1, fix="start") # TSS locations
  
  observeEvent(input$run_annotation, {
    tryCatch({
      req(consensus_peaks_rv())
      
      # Find nearest TSS for each peak
      peaks <- consensus_peaks_rv()
      nearest_idx <- nearest(peaks, tss_gr)
      peak_to_tss <- tss_gr[nearest_idx]
      
      # Calculate distance to TSS
      dist_to_tss <- start(peaks) - start(peak_to_tss)
      
      # Annotate with Entrez ID and Symbol
      entrez_ids <- names(peak_to_tss)
      gene_symbols <- mapIds(
        org.Hs.eg.db, 
        keys=entrez_ids, 
        column="SYMBOL", 
        keytype="ENTREZID", 
        multiVals="first"
      )
      
      # Construct annotation table
      ann_df <- data.frame(
        seqnames = as.character(seqnames(peaks)),
        start = start(peaks),
        end = end(peaks),
        width = width(peaks),
        strand = as.character(strand(peaks)),
        nearest_gene_id = entrez_ids,
        gene_symbol = gene_symbols,
        distance_to_tss = dist_to_tss
      )
      
      annotation_data(ann_df)
      showStatus("✅ Annotation complete")
    }, error=function(e){ 
      log_error(e, "run_annotation"); 
      showNotification("❌ Annotation failed", type="error") 
    })
  })
  
  output$peak_table <- renderDT({
    req(annotation_data())
    datatable(as.data.frame(annotation_data()))
  })
  
  output$download_peak_table <- downloadHandler(
    filename = "peak_annotations.csv",
    content = function(f) write.csv(as.data.frame(annotation_data()), f, row.names=FALSE)
  )
  
  output$pie_plot <- renderPlot({
    req(annotation_data())
    # Instead of plotAnnoPie (from ChIPseeker), we can make a custom plot, e.g. pie chart of gene_symbol counts
    ann_df <- annotation_data()
    gene_counts <- sort(table(ann_df$gene_symbol), decreasing=TRUE)
    pie(gene_counts[gene_counts > 0], main="Annotated Peaks (by nearest gene)")
  })
  


  
  # --- DAA (stub)
  observeEvent(input$run_daa_real, {
    tryCatch({
      req(input$count_matrix_file, input$metadata_file)
      
      counts <- read.csv(input$count_matrix_file$datapath, row.names = 1)
      metadata <- read.csv(input$metadata_file$datapath, row.names = 1)
      metadata <- metadata[colnames(counts), , drop = FALSE]
      
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts),
                                            colData = metadata,
                                            design = ~ condition)
      dds <- DESeq2::DESeq(dds)
      res <- DESeq2::results(dds)
      res_df <- as.data.frame(res[order(res$padj), ])
      res_df$peak <- rownames(res_df)
      
      uploaded_counts_rv(counts)
      uploaded_daa_rv(res_df)
      
      output$uploaded_daa_table <- renderDT({
        datatable(res_df, options = list(scrollX = TRUE))
      })
      
      output$download_uploaded_daa <- downloadHandler(
        filename = "daa_results.csv",
        content = function(file) write.csv(res_df, file, row.names = FALSE)
      )
      
      showStatus("✅ DAA complete.")
      
    }, error = function(e) {
      log_error(e, "run_daa_real")
      showNotification(paste("❌ DAA error:", conditionMessage(e)), type = "error")
    })
  })
  
  
  # --- PCA & UMAP (stub)
  observeEvent(input$run_pca, {
    tryCatch({
      req(uploaded_counts_rv(), input$metadata_file)
      
      counts <- uploaded_counts_rv()
      metadata <- read.csv(input$metadata_file$datapath, row.names = 1)
      log_counts <- log2(counts + 1)
      
      pca <- prcomp(t(log_counts), scale. = TRUE)
      pca_df <- as.data.frame(pca$x[, 1:2])
      pca_df$condition <- metadata[rownames(pca_df), "condition"]
      
      output$pca_plot <- renderPlotly({
        plot_ly(pca_df, x = ~PC1, y = ~PC2, color = ~condition, type = 'scatter', mode = 'markers') %>%
          layout(title = "PCA: PC1 vs PC2")
      })
      
      pca_data_rv(pca_df)
      showStatus("✅ PCA complete.")
      
    }, error = function(e) {
      log_error(e, "run_pca")
      showNotification(paste("❌ PCA error:", conditionMessage(e)), type = "error")
    })
  })
  observeEvent(input$run_umap, {
    tryCatch({
      req(uploaded_counts_rv(), input$metadata_file)
      
      counts <- uploaded_counts_rv()
      metadata <- read.csv(input$metadata_file$datapath, row.names = 1)
      log_counts <- log2(counts + 1)
      
      umap_res <- umap::umap(t(log_counts), n_neighbors = 5, min_dist = 0.3)
      umap_df <- as.data.frame(umap_res$layout)
      colnames(umap_df) <- c("UMAP1", "UMAP2")
      umap_df$condition <- metadata[rownames(umap_df), "condition"]
      
      output$umap_plot <- renderPlotly({
        plot_ly(umap_df, x = ~UMAP1, y = ~UMAP2, color = ~condition, type = 'scatter', mode = 'markers') %>%
          layout(title = "UMAP Plot")
      })
      
      umap_data_rv(umap_df)
      showStatus("✅ UMAP complete.")
      
    }, error = function(e) {
      log_error(e, "run_umap")
      showNotification(paste("❌ UMAP error:", conditionMessage(e)), type = "error")
    })
  })
  
  
  # --- Motif Enrichment (stub)
  observeEvent(input$run_motif, {
    tryCatch({
      req(consensus_peaks_rv())
      
      showStatus("📥 Fetching motifs from JASPAR2020")
      opts <- if (input$jaspar_family == "ALL") list() else list(species = 9606, family = input$jaspar_family)
      pwm_list <- TFBSTools::getMatrixSet(JASPAR2020, opts)
      
      peak_gr <- consensus_peaks_rv()
      showStatus(paste0("📏 Consensus peak count: ", length(peak_gr)))
      
      showStatus("🔍 Matching motifs to peaks...")
      matches <- motifmatchr::matchMotifs(pwm_list, peak_gr, genome = BSgenome.Hsapiens.UCSC.hg38)
      motif_mat <- motifmatchr::motifMatches(matches)
      showStatus("✅ Motif matching complete")
      
      showStatus(paste0("🧮 Matrix dimensions: ", nrow(motif_mat), " x ", ncol(motif_mat)))
      
      # Use colSums instead of rowSums to evaluate per-motif presence
      motif_scores <- Matrix::colSums(motif_mat)
      
      if (length(motif_scores) == 0 || all(motif_scores == 0)) {
        showNotification("⚠️ No motif enrichment detected (all zeros). Try more peaks or different TF family.", type = "warning")
        return()
      }
      
      motif_df <- data.frame(
        motif = names(motif_scores),
        count = motif_scores,
        perc_peaks = motif_scores / nrow(motif_mat) * 100
      )
      motif_df <- motif_df[order(-motif_df$count), ][1:min(20, nrow(motif_df)), ]
      
      motif_enrich_rv(motif_df)
      
      output$motif_enrich_plot <- renderPlotly({
        plot_ly(motif_df, x = ~motif, y = ~count, type = "bar") %>%
          layout(title = "Top Motifs", xaxis = list(title = ""), yaxis = list(title = "Peak Count"))
      })
      
      output$motif_enrich_table <- renderDT({
        datatable(motif_df)
      })
      
      output$download_motif <- downloadHandler(
        filename = function() paste0("motif_enrichment_", Sys.Date(), ".csv"),
        content = function(file) {
          write.csv(motif_df, file, row.names = FALSE)
        }
      )
      
      showStatus("✅ Motif enrichment complete.")
      
    }, error = function(e) {
      log_error(e, "run_motif")
      showNotification(paste("❌ Motif error:", conditionMessage(e)), type = "error")
    })
  })
  
  
  # --- Random Forest + ROC
  observeEvent(input$run_rf,{
    tryCatch({
      req(uploaded_counts_rv(), input$metadata_file)
      counts <- t(uploaded_counts_rv())
      meta <- read.csv(input$metadata_file$datapath, row.names=1)
      labels <- factor(meta[rownames(counts),"condition"])
      idx <- sample(seq_along(labels), length(labels)%/%2)
      rf <- randomForest(counts[idx,], labels[idx], ntree=500, importance=TRUE)
      preds <- predict(rf, counts[-idx,], type="prob")[,2]
      roc_obj <- roc(labels[-idx], preds)
      rf_results_rv(list(auc=auc(roc_obj), acc=mean(predict(rf, counts[-idx,])==labels[-idx])))
      rf_importance_rv(importance(rf))
      
      output$roc_plot <- renderPlotly({
        plot_ly(x=roc_obj$specificities, y=roc_obj$sensitivities, type="scatter", mode="lines") %>%
          layout(title=paste0("ROC (AUC=", round(auc(roc_obj),3),")"), xaxis=list(title="1-Specificity"), yaxis=list(title="Sensitivity"))
      })
      
      output$rf_metrics <- renderDT({
        df <- data.frame(Metric=c("AUC","Accuracy"), Value=c(rf_results_rv()$auc, rf_results_rv()$acc))
        datatable(df, options=list(dom='t'))
      })
      output$rf_varimp_plot <- renderPlotly({
        imp <- rf_importance_rv()[, "MeanDecreaseGini"]
        imp <- sort(imp, decreasing=TRUE)[1:20]
        plot_ly(x=names(imp), y=~imp, type="bar") %>% layout(xaxis=list(tickangle=-45))
      })
      
      showStatus("✅ RF complete")
    }, error=function(e){ log_error(e,"run_rf"); showNotification("❌ RF error",type="error") })
  })
  
  # --- Power Analysis
  observeEvent(input$run_power,{
    tryCatch({
      req(uploaded_counts_rv(), input$metadata_file)
      counts <- uploaded_counts_rv()
      logc <- log2(counts+1); reps <- input$power_replicates
      d <- mean(apply(logc,1,var), na.rm=TRUE)
      a <- input$power_alpha; eff <- input$power_effect_size
      df <- data.frame(Replicates=2:20)
      df$Power <- pnorm(sqrt(df$Replicates)*eff/sqrt(d)-qnorm(1-a/2)) + pnorm(-sqrt(df$Replicates)*eff/sqrt(d)-qnorm(1-a/2))
      df$CI_low <- pmax(0, df$Power - 1.96*sqrt(df$Power*(1-df$Power)/reps))
      df$CI_high <- pmin(1, df$Power + 1.96*sqrt(df$Power*(1-df$Power)/reps))
      power_df_rv(df)
      
      output$power_plot <- renderPlotly({
        plot_ly(df, x=~Replicates, y=~Power, type="scatter", mode="lines+markers") %>%
          layout(yaxis=list(range=c(0,1)))
      })
      output$power_table <- renderDT({ datatable(df) })
      output$download_power_table <- downloadHandler(
        filename = "power_table.csv",
        content = function(f) write.csv(power_df_rv(), f, row.names=FALSE)
      )
      showStatus("✅ Power analysis complete")
    }, error=function(e){ log_error(e,"run_power"); showNotification("❌ Power error",type="error") })
  })
  observeEvent(input$count_matrix_file, {
    tryCatch({
      req(input$count_matrix_file)
      counts <- read.csv(input$count_matrix_file$datapath, row.names = 1)
      uploaded_counts_rv(counts)
      showNotification("✅ Count matrix uploaded and ready.", type = "message")
    }, error = function(e) {
      log_error(e, "upload_counts")
      showNotification(paste("❌ Failed to load counts:", conditionMessage(e)), type = "error")
    })
  })
  observeEvent(input$metadata_file, {
    tryCatch({
      req(input$metadata_file)
      meta <- read.csv(input$metadata_file$datapath, row.names = 1)
      showNotification("✅ Metadata uploaded.", type = "message")
    }, error = function(e) {
      log_error(e, "upload_metadata")
      showNotification(paste("❌ Failed to load metadata:", conditionMessage(e)), type = "error")
    })
  })
  
  observeEvent(input$run_heatmap, {
    tryCatch({
      req(uploaded_counts_rv(), input$metadata_file)
      showNotification("📥 Retrieving uploaded count matrix...", type = "message")
      
      counts <- uploaded_counts_rv()
      log_counts <- log2(counts + 1)
      showNotification("🔍 Log-transformed counts", type = "message")
      
      # Variance filter
      row_vars <- apply(log_counts, 1, var)
      showNotification(paste("🔬 Variance calculated across", nrow(log_counts), "rows"), type = "message")
      
      top_var <- order(row_vars, decreasing = TRUE)[1:min(500, nrow(log_counts))]
      mat <- log_counts[top_var, , drop = FALSE]
      showNotification(paste("📊 Selected top", nrow(mat), "variable peaks"), type = "message")
      
      # Debug matrix size
      print(dim(mat))
      if (nrow(mat) < 2 || ncol(mat) < 2) {
        showNotification("⚠️ Not enough rows or columns to cluster", type = "error")
        return()
      }
      
      # Clustering
      showNotification("🔗 Performing clustering...", type = "message")
      row_dist <- dist(mat)
      row_clust <- hclust(row_dist)
      col_dist <- dist(t(mat))
      col_clust <- hclust(col_dist)
      showNotification("✅ Clustering complete", type = "message")
      
      # Reorder matrix
      heatmap_data <- as.matrix(mat[row_clust$order, col_clust$order])
      print("✅ Matrix reordered for heatmap")
      
      # Plotly heatmap
      output$heatmap_plot <- renderPlotly({
        plot_ly(
          z = heatmap_data,
          type = "heatmap",
          colors = colorRamp(c("navy", "white", "firebrick")),
          x = colnames(heatmap_data),
          y = rownames(heatmap_data)
        ) %>% layout(
          title = "Clustered Heatmap of Top Variable Peaks",
          xaxis = list(title = "Samples"),
          yaxis = list(title = "Peaks")
        )
      })
      
      showStatus("✅ Heatmap generated.")
      
    }, error = function(e) {
      msg <- paste("❌ Heatmap error:", conditionMessage(e))
      log_error(e, "run_heatmap")
      showNotification(msg, type = "error")
      print(msg)
    })
  })
  

}  # ← closes the server function

  

  


# ==== App Launch ====
shinyApp(ui, server)
