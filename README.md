# 🔬 JCAP_ATAC_SEQ APP

An interactive R Shiny app for analyzing and visualizing ATAC-seq peak data.  
Upload BED files, generate consensus peaks, annotate regions, run motif enrichment, perform machine learning, simulate power, and more — **no coding required**.

---

## 🚀 Features

- 📂 **Upload** ZIPs of BED files to generate consensus peaks
- 🧬 **Peak annotation** (nearest gene, distance to TSS, etc)
- 📊 **Motif enrichment** (using motifmatchr + JASPAR2020)
- 🧪 **Differential accessibility** (DAA) with DESeq2
- 📈 **PCA, UMAP, heatmaps** for exploratory analysis
- 🤖 **Random Forest** classifier (AUC, feature importance)
- 🔬 **Power analysis** for experimental design
- 💾 **Downloadable** results for all major outputs
- 🎨 **Anime-inspired UI** (with fairy_tail.css)
- 🖥️ Runs in RStudio, Shiny Server, or on HPC (Singularity)

---

## 🧪 Sample Data

To test the app’s full workflow, use the provided files:

| File                       | Purpose                         |
|----------------------------|---------------------------------|
| `focused_promoter_peaks.zip` | Consensus peaks (ZIP of BED)    |
| `simulated_counts (1).csv`   | Simulated peak × sample counts  |
| `simulated_metadata (1).csv` | Sample annotations (incl. group)|

**How to use:**  
Upload each file via the corresponding input in the sidebar.

---

## 📦 Requirements

You need R (≥ 4.3.x) and the following R packages:

```r
install.packages(c(
  "shiny", "shinyjs", "plotly", "DT", "markdown",
  "randomForest", "pROC", "enrichR","uwot", "umap"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "GenomicRanges", "GenomicFeatures", "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "TFBSTools",
  "JASPAR2020", "BSgenome.Hsapiens.UCSC.hg38", "motifmatchr"
))


🧬 How to Run the App
Option 1: Local RStudio
Download/clone this repository.

Open RStudio, set your working directory to the project folder:

and  run  this in your console

setwd("/path/to/ATAC_APP")
library(shiny)
runApp(".")
The app will open in your browser (e.g. http://localhost:8787 if using Shiny Server).

HPC / Singularity Deployment
Option 2: HPC (with SLURM & Singularity)

1.Build and run from the bash command line:

singularity build atac-shiny.sif Singularity.def
2.Submit your job with SLURM:

sbatch run_atac_app.sh

Monitor the output and error logs as specified in the batch script.

3.Forward port 8787 to your local machine for browser access:

ssh -L 8787:localhost:8787 user@cluster

Open the app in your browser:
http://localhost:8787
Tip:

Edit run_atac_app.sh to adjust resource requirements or output locations for your cluster.

See example SLURM script (run_atac_app.sh) included in the repo.


🎛️ How to Use the App
Consensus Peaks:
Upload a ZIP of BED files and click Make Consensus Peaks.

Peak Annotation:
Click Run Peak Annotation.

Results: Peak Annotation Table, Pie Chart.
🧬 Enrichr Pathway Analysis

Run Enrichr

Run DAA first

Select a contrast (e.g. Treated_vs_Control)

Choose a gene set (DAA_ALL, DAA_UP, or DAA_DOWN)

Click Run Enrichr

The app will:

Map ATAC peaks → genes

Send the gene set to Enrichr

Retrieve enriched GO, Reactome, KEGG, and regulatory pathways

Rank results by Combined Score (or Adjusted p-value when needed)

Filter out weak single-gene overlaps

Results appear in:

Enrichr Table (full pathway list)

Enrichr Barplot (top ranked biological programs)

This converts chromatin accessibility changes into interpretable biological pathways for downstream analysis and LLM triage.
Motif Enrichment:
Select a JASPAR family, click Run Motif Enrichment.

Results: Motif Table, Motif Barplot.

Counts & Metadata:
Upload count matrix (CSV) and metadata (CSV).
Click Run DAA on Uploaded Counts to see DESeq2 results.

.


Exploratory Visualizations:

PCA: Click Run PCA

UMAP: Click Run UMAP

Heatmap: Click Plot Heatmap

Random Forest Classifier:
Click Run Random Forest Classifier.

Results: AUC, Accuracy, Feature Importance.

Power Analysis:
Set effect size, FDR, and replicates. Click Run Power Analysis.

Downloadable Results:
All major tables have Download buttons.

| Tab Name                            | Description                                                                            |
| ----------------------------------- | -------------------------------------------------------------------------------------- |
| **README**                          | Usage guide, pipeline overview, and interpretation notes                               |
| **Peak Annotation Table**           | Annotated consensus peaks with nearest gene, TSS distance, and peak IDs (downloadable) |
| **Annotation Pie Chart**            | Visual distribution of peaks by nearest gene (top hits)                                |
| **Motif Enrichment Table**          | Transcription factor motifs matched to accessible regions                              |
| **Motif Enrichment Plot**           | Barplot of top enriched transcription factor motifs                                    |
| **DAA Results (per condition)**     | Differential accessibility results from DESeq2 for each contrast                       |
| **DAA Gene Sets (ALL / UP / DOWN)** | Enrichr-ready gene sets built from significant ATAC peaks                              |
| **Enrichr Pathways**                | Pathway enrichment of ATAC-derived gene sets (GO, Reactome, etc.)                      |
| **Enrichr Barplot**                 | Ranked pathway barplot (Combined Score / Adjusted p-value)                             |
| **PCA Plot**                        | Principal component analysis of ATAC peak counts                                       |
| **UMAP Plot**                       | UMAP projection of samples based on chromatin accessibility                            |
| **Heatmap**                         | Clustered heatmap of top variable peaks (log-scaled counts)                            |
| **RF Metrics**                      | Random forest classifier performance (AUC, accuracy)                                   |
| **Feature Importance**              | Peaks ranked by Gini importance from the RF model                                      |
| **ROC Curve**                       | Receiver-operating characteristic curve for classification                             |
| **Power Plot**                      | Statistical power vs. number of replicates                                             |
| **Power Table**                     | Power estimates with confidence intervals (downloadable)                               |


🛠️ Developer Notes
Error logging: All errors are appended to error_log.txt.

Logs: Use email_log.R to manage or email logs if desired.
ATAC_APP/
├── app.R                    # Main Shiny app
├── run.sh                   # Launch script (local)
├── run_atac_app.sh          # SLURM batch script for HPC
├── email_log.R              # Log emailer
├── error_log.txt            # Error logs
├── www/
│   └── fairy_tail.css       # Themed UI
├── sample_data/
│   └── focused_promoter_peaks.zip
│   └── simulated_counts (1).csv
│   └── simulated_metadata (1).csv
├── logs/                    # Archived logs
├── Singularity.def          # Singularity definition
└── README.md


🧠 FAQ
Can I use BAM files?
No — only BED/peak-level data are supported. Use MACS2 or similar to call peaks first.

Is this cloud deployable?
Not trivially — Bioconductor packages and large genome data require persistent local storage. HPC, local, or Docker/Singularity are recommended.

How big can my datasets be?
The app is robust for datasets up to ~10k peaks × 100 samples. Larger sets may need more RAM.
