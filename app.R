# app.R (multi-species annotation + motif update)
library(shiny); library(shinyjs); library(plotly); library(DT)

# Core genolibrary(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)

library(GenomicRanges)
library(GenomicFeatures)
library(uwot); library(randomForest); library(pROC)
library(umap)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(BSgenome.Drerio.UCSC.danRer11)

library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
library(BSgenome.Dmelanogaster.UCSC.dm6)

# Motifs
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)

# We will load species-specific packages lazily
# (TxDb.* , org.*.eg.db, BSgenome.*) based on user selection

options(shiny.maxRequestSize = 200*1024^2)

# ---- Multi-species helpers (explicit, robust) ----
get_txdb <- function(sp) {
  switch(sp,
         human     = TxDb.Hsapiens.UCSC.hg38.knownGene,
         mouse     = TxDb.Mmusculus.UCSC.mm10.knownGene,
         zebrafish = TxDb.Drerio.UCSC.danRer11.refGene,
         fly       = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
         TxDb.Hsapiens.UCSC.hg38.knownGene
  )
}

get_orgdb <- function(sp) {
  switch(sp,
         human     = org.Hs.eg.db,
         mouse     = org.Mm.eg.db,
         zebrafish = org.Dr.eg.db,
         fly       = org.Dm.eg.db,
         org.Hs.eg.db
  )
}

get_bsgenome <- function(sp) {
  switch(sp,
         human     = BSgenome.Hsapiens.UCSC.hg38,
         mouse     = BSgenome.Mmusculus.UCSC.mm10,
         zebrafish = BSgenome.Drerio.UCSC.danRer11,
         fly       = BSgenome.Dmelanogaster.UCSC.dm6,
         BSgenome.Hsapiens.UCSC.hg38
  )
}


 

get_jaspar_taxid <- function(sp) {
  switch(sp,
         human = 9606L,
         mouse = 10090L,
         zebrafish = 7955L,
         fly = 7227L,
         9606L
  )
}

# Build (TxDb, genes, TSS, orgdb, bsgenome, jaspar_taxid) bundle


# ---------- UI ----------
ui <- fluidPage(
  includeCSS("www/fairy_tail.css"),
  useShinyjs(),
  titlePanel("JCAP_ATAC_SEQ APP"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "species", "Species",
        choices = c(
          "Human (hg38)" = "human",
          "Mouse (mm10)" = "mouse",
          "Zebrafish (danRer11)" = "zebrafish",
          "Fly (dm6)" = "fly"
        ),
        selected = "human"
      ),
      fileInput("zipfile", "Upload ZIP of BED files for consensus peaks", accept = ".zip"),
      actionButton("make_peaks", "Make Consensus Peaks"),
      actionButton("run_annotation", "Run Peak Annotation"),
      actionButton("run_motif", "Run Motif Enrichment"),
      selectInput("jaspar_family", "TF Family (JASPAR)", choices = c("ALL","MAF","TP53","C2H2")),
      fileInput("count_matrix_file", "Upload Count Matrix CSV", accept = ".csv"),
      fileInput("metadata_file", "Upload Metadata CSV", accept = ".csv"),
      actionButton("run_daa_real", "Run DAA on Uploaded Counts"),
      actionButton("run_enrichr_btn", "Run Enrichr (Selected Contrast)"),
      selectInput(
        "enrichr_db",
        "Enrichr database",
        choices = c(
          "Reactome_2022",
          "GO_Biological_Process_2021",
          "KEGG_2021_Human",
          "WikiPathways_2023_Human",
          "MSigDB_Hallmark_2020"
        ),
        selected = "Reactome_2022"
      ),
      
      
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
        tabPanel(
          "DAA Results (per condition)",
          selectInput("contrast_pick", "Select contrast", choices = c("Run DAA first" = "")),
          hr(),
          h4("Full DAA table"),
          DTOutput("daa_table_selected"),
          downloadButton("dl_daa_selected", "Download full DAA CSV"),
          hr(),
          h4("Gene sets for Enrichment Triage App"),
          h5("DAA_ALL"),
          DTOutput("gs_all_selected"),
          downloadButton("dl_all_selected", "Download DAA_ALL geneset"),
          hr(),
          h5("DAA_UP"),
          DTOutput("gs_up_selected"),
          downloadButton("dl_up_selected", "Download DAA_UP geneset"),
          hr(),
          h5("DAA_DOWN"),
          DTOutput("gs_down_selected"),
          downloadButton("dl_down_selected", "Download DAA_DOWN geneset"),
          hr(),
          downloadButton("dl_bundle_selected", "Download ALL (ALL+UP+DOWN) as one CSV")
        ),
        tabPanel(
          "Enrichr (per contrast)",
          selectInput("enrichr_set_pick", "Gene set", choices = c("ALL", "UP", "DOWN"), selected = "ALL"),
          numericInput("enrichr_top_n", "Top terms to show", value = 15, min = 5, max = 50),
          
          tabsetPanel(
            tabPanel(
              "Barplot",
              plotlyOutput("enrichr_barplot", height = "650px")
            ),
            tabPanel(
              "Tables",
              h5("ALL"),
              DTOutput("enrichr_all_tbl"),
              downloadButton("dl_enrichr_all", "Download ALL"),
              hr(),
              h5("UP"),
              DTOutput("enrichr_up_tbl"),
              downloadButton("dl_enrichr_up", "Download UP"),
              hr(),
              h5("DOWN"),
              DTOutput("enrichr_down_tbl"),
              downloadButton("dl_enrichr_down", "Download DOWN"),
              hr(),
              downloadButton("dl_enrichr_bundle", "Download ALL+UP+DOWN bundle")
            )
          )
        ),
        
        
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

# ---------- SERVER ----------
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
  enrichr_by_contrast_rv <- reactiveVal(list())  # contrast -> list(db=..., all=..., up=..., down=...)
  
  daa_by_contrast_rv <- reactiveVal(list())   # contrast -> full DAA res_df
  gene_sets_rv       <- reactiveVal(list())   # contrast -> list(all/up/down geneset dfs)
  
  showStatus <- function(msg) showNotification(msg, type="message", duration=5)
  log_error <- function(e, ctx) write(paste(Sys.time(), ctx, conditionMessage(e)), file="error_log.txt", append=TRUE)
  
  # Build a species bundle reactive that updates when species changes
  
  # Convert "A;B;C" into a character vector safely
  split_genes <- function(s) {
    if (is.null(s) || is.na(s) || !nzchar(s)) return(character())
    g <- unlist(strsplit(s, ";", fixed = TRUE))
    g <- unique(trimws(g))
    g[nzchar(g)]
  }
  
  
  prep_enrichr_for_plot <- function(tbl, topn = 15, min_overlap_genes = 2) {
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) return(NULL)

  # Detect columns (Enrichr naming varies a bit)
  term_col <- if ("Term" %in% colnames(tbl)) "Term" else colnames(tbl)[1]

  # Prefer adjusted p if present
  p_col <- if ("Adjusted.P.value" %in% colnames(tbl)) "Adjusted.P.value"
  else if ("Adjusted P-value" %in% colnames(tbl)) "Adjusted P-value"
  else if ("P.value" %in% colnames(tbl)) "P.value"
  else if ("P-value" %in% colnames(tbl)) "P-value"
  else NULL

  # Combined score if present
  cs_col <- if ("Combined.Score" %in% colnames(tbl)) "Combined.Score"
  else if ("Combined Score" %in% colnames(tbl)) "Combined Score"
  else NULL

  # Overlap column is usually like "3/250"
  ov_col <- if ("Overlap" %in% colnames(tbl)) "Overlap" else NULL

  # Must have p-values to plot
  if (is.null(p_col)) return(NULL)

  df <- tbl

  # Parse overlap numerator for filtering
  if (!is.null(ov_col)) {
    ov_num <- suppressWarnings(as.integer(sub("/.*$", "", df[[ov_col]])))
    ov_num[is.na(ov_num)] <- 0L
    df$OverlapN <- ov_num
    df <- df[df$OverlapN >= min_overlap_genes, , drop = FALSE]
  }

  if (nrow(df) == 0) return(NULL)

  # Ensure p-values sane; avoid -Inf
  df[[p_col]] <- suppressWarnings(as.numeric(df[[p_col]]))
  df[[p_col]][is.na(df[[p_col]])] <- 1
  df$mlog10p <- -log10(pmax(df[[p_col]], 1e-300))

  # Ranking:
  # 1) If Combined Score exists, rank by it (desc), tie-break by p (asc)
  # 2) else rank by p (asc)
  if (!is.null(cs_col)) {
    df[[cs_col]] <- suppressWarnings(as.numeric(df[[cs_col]]))
    df[[cs_col]][is.na(df[[cs_col]])] <- -Inf
    ord <- order(-df[[cs_col]], df[[p_col]])
  } else {
    ord <- order(df[[p_col]])
  }
  df <- df[ord, , drop = FALSE]

  # Take top N
  df <- df[1:min(topn, nrow(df)), , drop = FALSE]

  list(
    df = df,
    term_col = term_col,
    p_col = p_col,
    cs_col = cs_col,
    ov_col = ov_col
  )
}

  # Run Enrichr safely; returns a DF (or empty DF)
  run_enrichr_safe <- function(genes, db) {
    genes <- unique(trimws(genes))
    genes <- genes[nzchar(genes)]
    if (length(genes) < 5) return(data.frame())  # too few genes
    
    res <- enrichR::enrichr(genes, db)
    out <- res[[db]]
    if (is.null(out)) return(data.frame())
    out
  }
  
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
      
      cons <- reduce(do.call(c, seqs))
      names(cons) <- paste0("peak_", seq_along(cons))
      consensus_peaks_rv(cons)
      
      showStatus(paste("✅ Consensus peaks:", length(consensus_peaks_rv())))
    }, error = function(e) {
      log_error(e, "make_peaks")
      showNotification("❌ Consensus peak generation failed", type = "error")
    })
  })
  
  
  observeEvent(input$run_annotation, {
    tryCatch({
      req(consensus_peaks_rv())
      
      sp <- switch(input$species,
                   human = "human", mouse = "mouse",
                   zebrafish = "zebrafish", fly = "fly", "human")
      
      txdb  <- get_txdb(sp)
      orgdb <- get_orgdb(sp)
      
      peaks <- consensus_peaks_rv()
      
      genes_gr <- GenomicFeatures::genes(txdb)
      
      # make TSS (keep strand so "start" is strand-aware)
      tss_gr <- promoters(genes_gr, upstream = 0, downstream = 1)
      
      nearest_idx <- nearest(peaks, tss_gr, ignore.strand = TRUE)
      
      # initialize outputs safely
      entrez_or_geneid <- rep(NA_character_, length(peaks))
      dist_to_tss      <- rep(NA_integer_,  length(peaks))
      
      ok <- !is.na(nearest_idx)
      if (any(ok)) {
        nearest_tss <- tss_gr[nearest_idx[ok]]
        dist_to_tss[ok] <- start(peaks[ok]) - start(nearest_tss)
        
        
        # Robust: TxDb gene IDs usually live in mcols(..)$gene_id
        if ("gene_id" %in% colnames(mcols(nearest_tss))) {
          entrez_or_geneid[ok] <- as.character(mcols(nearest_tss)$gene_id)
        } else {
          entrez_or_geneid[ok] <- as.character(names(nearest_tss))
        } }
        
      
      # ---- auto-detect keytype for mapIds based on ID patterns ----
      detect_keytype <- function(ids) {
        ids <- ids[!is.na(ids)]
        if (!length(ids)) return("ENTREZID")
        
        if (any(grepl("^FBgn", ids))) return("FLYBASE")
        if (any(grepl("^ENS",  ids))) return("ENSEMBL")
        if (all(grepl("^[0-9]+$", ids))) return("ENTREZID")
        
        # fallback
        "ENTREZID"
      }
      
      kt <- detect_keytype(entrez_or_geneid)
      
      gene_symbols <- rep(NA_character_, length(peaks))
      if (any(ok)) {
        gene_symbols[ok] <- unname(AnnotationDbi::mapIds(
          orgdb,
          keys      = entrez_or_geneid[ok],
          column    = "SYMBOL",
          keytype   = kt,
          multiVals = "first"
        ))
      }
      
      ann_df <- data.frame(
        seqnames        = as.character(seqnames(peaks)),
        start           = start(peaks),
        end             = end(peaks),
        width           = width(peaks),
        strand          = as.character(strand(peaks)),
        nearest_gene_id = entrez_or_geneid,
        gene_symbol     = gene_symbols,
        distance_to_tss = dist_to_tss,
        stringsAsFactors = FALSE
      )
      
      ann_df$peak_id <- names(peaks)
      
      
      annotation_data(ann_df)
      showNotification(paste0("✅ Annotation complete (keytype=", kt, ")"), type="message")
      
    }, error = function(e) {
      write(paste(Sys.time(), "run_annotation", conditionMessage(e)),
            file = "error_log.txt", append = TRUE)
      showNotification(paste("❌ Annotation failed:", conditionMessage(e)), type = "error")
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
    ann_df <- annotation_data()
    gene_counts <- sort(table(ann_df$gene_symbol), decreasing=TRUE)
    # show top 20 to keep readable
    if (length(gene_counts) > 20) gene_counts <- gene_counts[1:20]
    pie(gene_counts, main="Annotated Peaks (Top Genes by Nearest TSS)")
  })
  
  # --- DAA (ATAC-safe, per-contrast, with robust geneset export) ---
  observeEvent(input$run_daa_real, {
    tryCatch({
      req(input$count_matrix_file, input$metadata_file)
      req(annotation_data())  # need peak_id -> gene mapping
      
      counts   <- read.csv(input$count_matrix_file$datapath, row.names = 1, check.names = FALSE)
      metadata <- read.csv(input$metadata_file$datapath, row.names = 1, check.names = FALSE)
      
      # ---- HARD sample alignment checks ----
      if (!all(colnames(counts) %in% rownames(metadata))) {
        missing <- setdiff(colnames(counts), rownames(metadata))
        stop("Metadata is missing these samples from the count matrix: ", paste(missing, collapse = ", "))
      }
      metadata <- metadata[colnames(counts), , drop = FALSE]
      if (!identical(rownames(metadata), colnames(counts))) {
        stop("Sample order mismatch: rownames(metadata) must exactly match colnames(counts).")
      }
      
      # ---- Ensure + clean condition ----
      req("condition" %in% colnames(metadata))
      metadata$condition <- factor(trimws(as.character(metadata$condition)))
      if (nlevels(metadata$condition) < 2) {
        stop("Need at least 2 condition levels. Found: ", paste(levels(metadata$condition), collapse = ", "))
      }
      
      # Prefer Control as reference if present
      if ("Control" %in% levels(metadata$condition)) {
        metadata$condition <- relevel(metadata$condition, ref = "Control")
      }
      
      # ---- ATAC peak filtering (CRITICAL) ----
      keep <- rowSums(counts >= 10) >= 2
      counts <- counts[keep, , drop = FALSE]
      if (nrow(counts) < 10) {
        stop("After filtering, too few peaks remain (", nrow(counts), "). Lower the filter threshold for this dataset.")
      }
      
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(counts),
        colData   = metadata,
        design    = ~ condition
      )
      
      # ---- Freeze size factors (prevents DESeq2 from normalizing away global shifts) ----
      DESeq2::sizeFactors(dds) <- rep(1, ncol(dds))
      
      # Run DESeq without estimating size factors (they're already set)
      dds <- DESeq2::DESeq(dds, sfType = "poscounts")
      
      
      # Baseline + contrasts
      ref_level  <- levels(metadata$condition)[1]
      other_lvls <- setdiff(levels(metadata$condition), ref_level)
      
      # ---- Annotation mapping ----
      ann <- annotation_data()
      if (!"peak_id" %in% colnames(ann)) {
        stop("annotation_data() is missing peak_id. Re-run annotation after consensus peaks are named peak_#.")
      }
      
      # helper: enrichr-style gene set table
      make_geneset_row <- function(term, peaks, padj_value = NA_real_, score = NA_real_) {
        genes <- unique(na.omit(ann$gene_symbol[match(peaks, ann$peak_id)]))
        genes <- genes[genes != ""]
        data.frame(
          Term = term,
          `Adjusted P-value` = padj_value,
          `Combined Score`   = score,
          Genes             = if (length(genes)) paste(genes, collapse = ";") else "",
          stringsAsFactors  = FALSE
        )
      }
      
      # Geneset fallback size if no padj hits
      TOP_N_FALLBACK <- 500
      
      daa_list <- list()
      gs_list  <- list()
      
      for (lvl in other_lvls) {
        contrast_name <- paste0(lvl, "_vs_", ref_level)
        
        res <- DESeq2::results(dds, contrast = c("condition", lvl, ref_level))
        res_df <- as.data.frame(res)
        res_df$peak_id <- rownames(res_df)
        
        # sort by padj then pvalue
        res_df <- res_df[order(res_df$padj, res_df$pvalue), , drop = FALSE]
        daa_list[[contrast_name]] <- res_df
        
        res_df2 <- res_df[!is.na(res_df$padj), ]
        
        # Take top N peaks by significance
        TOP_N <- 1000   # 500–2000 is usually perfect
        res_top <- res_df2[order(res_df2$padj), ][1:min(TOP_N, nrow(res_df2)), ]
        
        sig <- res_top
        sig_up   <- subset(res_top, log2FoldChange > 0)
        sig_down <- subset(res_top, log2FoldChange < 0)
        
        
        # Fallback if nothing significant: take top N by pvalue (keeps your pipeline moving)
        if (nrow(sig) == 0) {
          sig <- res_df[order(res_df$pvalue), , drop = FALSE]
          sig <- sig[!is.na(sig$pvalue), , drop = FALSE]
          sig <- sig[1:min(TOP_N_FALLBACK, nrow(sig)), , drop = FALSE]
          showNotification(
            paste0("⚠️ ", contrast_name, ": no padj<0.05 peaks. Using top ", nrow(sig), " by p-value for genesets."),
            type = "warning"
          )
        }
        
        sig_up   <- subset(sig, log2FoldChange > 0)
        sig_down <- subset(sig, log2FoldChange < 0)
        
        minp_all  <- if (nrow(sig)      > 0) min(sig$padj, na.rm = TRUE) else NA_real_
        minp_up   <- if (nrow(sig_up)   > 0) min(sig_up$padj, na.rm = TRUE) else NA_real_
        minp_down <- if (nrow(sig_down) > 0) min(sig_down$padj, na.rm = TRUE) else NA_real_
        
        gs_all  <- make_geneset_row(paste0("DAA_ALL_",  contrast_name), sig$peak_id,      minp_all)
        gs_up   <- make_geneset_row(paste0("DAA_UP_",   contrast_name), sig_up$peak_id,   minp_up)
        gs_down <- make_geneset_row(paste0("DAA_DOWN_", contrast_name), sig_down$peak_id, minp_down)
        
        gs_list[[contrast_name]] <- list(all = gs_all, up = gs_up, down = gs_down)
        
        # quick status
        showNotification(
          paste0("✅ ", contrast_name,
                 " | peaks kept=", nrow(counts),
                 " | sig(padj<0.05)=", sum(!is.na(res_df$padj) & res_df$padj < 0.05)),
          type = "message"
        )
      }
      
      daa_by_contrast_rv(daa_list)
      gene_sets_rv(gs_list)
      
      showStatus(paste0("✅ DAA complete for ", length(other_lvls), " contrast(s)."))
      
    }, error = function(e) {
      log_error(e, "run_daa_real_multi")
      showNotification(paste("❌ DAA error:", conditionMessage(e)), type = "error")
    })
  })
  
  
  # --- PCA & UMAP (unchanged stubs) ---
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
  
  # --- Motif Enrichment (species-aware, using helpers) ---
  observeEvent(input$run_motif, {
    tryCatch({
      req(consensus_peaks_rv())
      
      # Map UI value to helper key (adjust if your selectInput uses different values)
      sp <- switch(input$species,
                   human = "human",
                   mouse = "mouse",
                   zebrafish = "zebrafish",
                   fly = "fly",
                   "human")
      
      # Get genome + JASPAR species id explicitly
      bsgenome <- get_bsgenome(sp)
      if (is.null(bsgenome)) {
        showNotification("⚠️ BSgenome for selected species is not installed; motif scanning cannot run.",
                         type = "warning")
        return()
      }
      taxid <- get_jaspar_taxid(sp)
      
      showStatus("📥 Fetching motifs from JASPAR2020")
      opts <- list(species = taxid)
      if (input$jaspar_family != "ALL") opts$family <- input$jaspar_family
      pwm_list <- TFBSTools::getMatrixSet(JASPAR2020, opts)
      
      peak_gr <- consensus_peaks_rv()
      showStatus(paste0("📏 Consensus peak count: ", length(peak_gr)))
      
      showStatus("🔍 Matching motifs to peaks...")
      matches   <- motifmatchr::matchMotifs(pwm_list, peak_gr, genome = bsgenome)
      motif_mat <- motifmatchr::motifMatches(matches)
      showStatus("✅ Motif matching complete")
      
      scores <- Matrix::colSums(motif_mat)
      
      if (length(scores) == 0 || all(scores == 0)) {
        showNotification("⚠️ No motif enrichment detected (all zeros). Try more peaks or a different TF family.",
                         type = "warning")
        return()
      }
      
      motif_df <- data.frame(
        motif      = names(scores),
        count      = as.numeric(scores),
        perc_peaks = as.numeric(scores) / nrow(motif_mat) * 100,
        stringsAsFactors = FALSE
      )
      motif_df <- motif_df[order(-motif_df$count), ][1:min(20, nrow(motif_df)), ]
      
      motif_enrich_rv(motif_df)
      
      output$motif_enrich_plot <- renderPlotly({
        plot_ly(motif_df, x = ~motif, y = ~count, type = "bar") %>%
          layout(title = "Top Motifs", xaxis = list(title = ""),
                 yaxis = list(title = "Peak Count"))
      })
      
      output$motif_enrich_table <- renderDT({
        DT::datatable(motif_df)
      })
      
      output$download_motif <- downloadHandler(
        filename = function() paste0("motif_enrichment_", Sys.Date(), ".csv"),
        content  = function(file) write.csv(motif_df, file, row.names = FALSE)
      )
      
      showStatus("✅ Motif enrichment complete.")
    }, error = function(e) {
      log_error(e, "run_motif")
      showNotification(paste("❌ Motif error:", conditionMessage(e)), type = "error")
    })
  })
  
  
  # --- Random Forest + ROC (unchanged) ---
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
  
  # --- Power Analysis (unchanged) ---
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
  
  # Upload hooks (unchanged)
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
  
  # Heatmap (unchanged)
  observeEvent(input$run_heatmap, {
    tryCatch({
      req(uploaded_counts_rv(), input$metadata_file)
      showNotification("📥 Retrieving uploaded count matrix...", type = "message")
      
      counts <- uploaded_counts_rv()
      log_counts <- log2(counts + 1)
      showNotification("🔍 Log-transformed counts", type = "message")
      
      row_vars <- apply(log_counts, 1, var)
      showNotification(paste("🔬 Variance calculated across", nrow(log_counts), "rows"), type = "message")
      
      top_var <- order(row_vars, decreasing = TRUE)[1:min(500, nrow(log_counts))]
      mat <- log_counts[top_var, , drop = FALSE]
      showNotification(paste("📊 Selected top", nrow(mat), "variable peaks"), type = "message")
      
      if (nrow(mat) < 2 || ncol(mat) < 2) {
        showNotification("⚠️ Not enough rows or columns to cluster", type = "error")
        return()
      }
      
      showNotification("🔗 Performing clustering...", type = "message")
      row_dist <- dist(mat)
      row_clust <- hclust(row_dist)
      col_dist <- dist(t(mat))
      col_clust <- hclust(col_dist)
      showNotification("✅ Clustering complete", type = "message")
      
      heatmap_data <- as.matrix(mat[row_clust$order, col_clust$order])
      
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
  
  
  
  observeEvent(daa_by_contrast_rv(), {
    contrasts <- names(daa_by_contrast_rv())
    if (length(contrasts) == 0) return()
    
    updateSelectInput(
      session,
      "contrast_pick",
      choices = setNames(contrasts, contrasts),
      selected = contrasts[1]
    )
  })
  selected_contrast <- reactive({
    req(input$contrast_pick)
    req(input$contrast_pick != "")
    input$contrast_pick
  })
  
  output$daa_table_selected <- renderDT({
    ct <- selected_contrast()
    datatable(daa_by_contrast_rv()[[ct]], options = list(scrollX = TRUE, pageLength = 25))
  })
  
  output$gs_all_selected <- renderDT({
    ct <- selected_contrast()
    datatable(gene_sets_rv()[[ct]]$all, options = list(dom = "t"))
  })
  
  output$gs_up_selected <- renderDT({
    ct <- selected_contrast()
    datatable(gene_sets_rv()[[ct]]$up, options = list(dom = "t"))
  })
  
  output$gs_down_selected <- renderDT({
    ct <- selected_contrast()
    datatable(gene_sets_rv()[[ct]]$down, options = list(dom = "t"))
  })
  
  output$dl_daa_selected <- downloadHandler(
    filename = function() paste0("DAA_", selected_contrast(), "_full.csv"),
    content = function(file) write.csv(daa_by_contrast_rv()[[selected_contrast()]], file, row.names = FALSE)
  )
  
  output$dl_all_selected <- downloadHandler(
    filename = function() paste0("DAA_ALL_", selected_contrast(), "_geneset.csv"),
    content = function(file) write.csv(gene_sets_rv()[[selected_contrast()]]$all, file, row.names = FALSE)
  )
  
  output$dl_up_selected <- downloadHandler(
    filename = function() paste0("DAA_UP_", selected_contrast(), "_geneset.csv"),
    content = function(file) write.csv(gene_sets_rv()[[selected_contrast()]]$up, file, row.names = FALSE)
  )
  
  output$dl_down_selected <- downloadHandler(
    filename = function() paste0("DAA_DOWN_", selected_contrast(), "_geneset.csv"),
    content = function(file) write.csv(gene_sets_rv()[[selected_contrast()]]$down, file, row.names = FALSE)
  )
  
  output$dl_bundle_selected <- downloadHandler(
    filename = function() paste0("DAA_bundle_", selected_contrast(), "_genesets.csv"),
    content = function(file) {
      gs <- gene_sets_rv()[[selected_contrast()]]
      write.csv(rbind(gs$all, gs$up, gs$down), file, row.names = FALSE)
    }
  )
  
  observeEvent(input$run_enrichr_btn, {
    tryCatch({
      ct <- selected_contrast()        # uses your reactive() that req()s contrast_pick
      req(gene_sets_rv()[[ct]])
      db <- input$enrichr_db
      
      gs <- gene_sets_rv()[[ct]]
      
      # Genes live in the Enrichr-style table row you created
      genes_all  <- split_genes(gs$all$Genes[1])
      genes_up   <- split_genes(gs$up$Genes[1])
      genes_down <- split_genes(gs$down$Genes[1])
      
      showStatus(paste0("🌐 Running Enrichr: ", ct, " | DB=", db))
      
      en_all  <- run_enrichr_safe(genes_all,  db)
      en_up   <- run_enrichr_safe(genes_up,   db)
      en_down <- run_enrichr_safe(genes_down, db)
      
      # store (keep previous contrasts)
      store <- enrichr_by_contrast_rv()
      store[[ct]] <- list(db = db, all = en_all, up = en_up, down = en_down)
      enrichr_by_contrast_rv(store)
      
      showStatus("✅ Enrichr complete.")
      
      # friendly warning if empty
      if (nrow(en_all) == 0 && nrow(en_up) == 0 && nrow(en_down) == 0) {
        showNotification("⚠️ Enrichr returned no results (often too few genes or no internet access).", type="warning")
      }
      
    }, error = function(e) {
      log_error(e, "run_enrichr_btn")
      showNotification(paste("❌ Enrichr error:", conditionMessage(e)), type="error")
    })
  })
  output$enrichr_all_tbl <- renderDT({
    ct <- selected_contrast()
    x <- enrichr_by_contrast_rv()[[ct]]$all
    req(!is.null(x))
    DT::datatable(x, options = list(scrollX = TRUE, pageLength = 25))
  })
  
  output$enrichr_up_tbl <- renderDT({
    ct <- selected_contrast()
    x <- enrichr_by_contrast_rv()[[ct]]$up
    req(!is.null(x))
    DT::datatable(x, options = list(scrollX = TRUE, pageLength = 25))
  })
  
  output$enrichr_down_tbl <- renderDT({
    ct <- selected_contrast()
    x <- enrichr_by_contrast_rv()[[ct]]$down
    req(!is.null(x))
    DT::datatable(x, options = list(scrollX = TRUE, pageLength = 25))
  })
  output$enrichr_barplot <- renderPlotly({
    ct <- selected_contrast()
    req(enrichr_by_contrast_rv()[[ct]])
    
    set_pick <- input$enrichr_set_pick
    topn <- input$enrichr_top_n
    
    tbl <- switch(set_pick,
                  "ALL"  = enrichr_by_contrast_rv()[[ct]]$all,
                  "UP"   = enrichr_by_contrast_rv()[[ct]]$up,
                  "DOWN" = enrichr_by_contrast_rv()[[ct]]$down)
    
    # Filter tiny overlaps + pick adjusted p + rank by combined score
    prepped <- prep_enrichr_for_plot(tbl, topn = topn, min_overlap_genes = 2)
    
    if (is.null(prepped)) {
      showNotification("⚠️ No Enrichr rows to plot after filtering (try lowering overlap filter or increasing gene set size).",
                       type="warning")
      return(plotly::plot_ly())
    }
    
    df <- prepped$df
    tcol <- prepped$term_col
    pcol <- prepped$p_col
    cscol <- prepped$cs_col
    ovcol <- prepped$ov_col
    
    # Make a nicer hover label
    hover_txt <- df[[tcol]]
    hover_txt <- paste0(
      "<b>", df[[tcol]], "</b>",
      "<br>-log10(p): ", round(df$mlog10p, 3),
      "<br>", if (!is.null(pcol)) paste0(pcol, ": ", signif(df[[pcol]], 4)) else "",
      if (!is.null(cscol)) paste0("<br>", cscol, ": ", round(df[[cscol]], 3)) else "",
      if (!is.null(ovcol)) paste0("<br>Overlap: ", df[[ovcol]]) else ""
    )
    
    plotly::plot_ly(
      df,
      x = ~mlog10p,
      y = ~reorder(df[[tcol]], mlog10p),
      type = "bar",
      orientation = "h",
      text = hover_txt,
      hoverinfo = "text"
    ) %>% plotly::layout(
      title = paste0("Enrichr ", set_pick, " | ", ct, " | ", enrichr_by_contrast_rv()[[ct]]$db,
                     " (rank=", if (!is.null(cscol)) "CombinedScore" else "p-value",
                     ", minOverlap≥2)"),
      xaxis = list(title = "-log10(Adjusted p or p)"),
      yaxis = list(title = "")
    )
  })
  
  output$dl_enrichr_all <- downloadHandler(
    filename = function() paste0("Enrichr_ALL_", selected_contrast(), "_", input$enrichr_db, ".csv"),
    content = function(file) {
      ct <- selected_contrast()
      write.csv(enrichr_by_contrast_rv()[[ct]]$all, file, row.names = FALSE)
    }
  )
  
  output$dl_enrichr_up <- downloadHandler(
    filename = function() paste0("Enrichr_UP_", selected_contrast(), "_", input$enrichr_db, ".csv"),
    content = function(file) {
      ct <- selected_contrast()
      write.csv(enrichr_by_contrast_rv()[[ct]]$up, file, row.names = FALSE)
    }
  )
  
  output$dl_enrichr_down <- downloadHandler(
    filename = function() paste0("Enrichr_DOWN_", selected_contrast(), "_", input$enrichr_db, ".csv"),
    content = function(file) {
      ct <- selected_contrast()
      write.csv(enrichr_by_contrast_rv()[[ct]]$down, file, row.names = FALSE)
    }
  )
  
  output$dl_enrichr_bundle <- downloadHandler(
    filename = function() paste0("Enrichr_BUNDLE_", selected_contrast(), "_", input$enrichr_db, ".csv"),
    content = function(file) {
      ct <- selected_contrast()
      x <- enrichr_by_contrast_rv()[[ct]]
      # add a label column so bundle is readable
      all  <- if (nrow(x$all))  cbind(Set="ALL",  x$all)  else data.frame(Set="ALL")
      up   <- if (nrow(x$up))   cbind(Set="UP",   x$up)   else data.frame(Set="UP")
      down <- if (nrow(x$down)) cbind(Set="DOWN", x$down) else data.frame(Set="DOWN")
      write.csv(rbind(all, up, down), file, row.names = FALSE)
    }
  )
  
  
}  # server

shinyApp(ui, server)
  


  


