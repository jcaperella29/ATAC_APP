

🔬 ATAC-seq Peak Annotation & Enrichment Viewer
An interactive R Shiny app for analyzing and visualizing ATAC-seq peak data. Upload BED files, generate consensus peaks, annotate regions, perform enrichment, train machine learning models, simulate power, and more — no coding required.

💡 Built in R with glowy CSS and rigorously tested on real ATAC-seq datasets.

📚 Table of Contents
📸 App Screenshots

🚀 Features

🧪 Sample Data

📦 Requirements

🧬 How to Run the App

🎛️ How to Use the App

📊 Tab Overview

🛠️ Developer Notes

📁 Folder Structure

👨‍🔬 Citation / Credit

🧠 FAQ

📸 App Screenshots
🧬 Peak Annotation Pie Chart


📑 Enrichment Table


📊 Peak Annotation Table


🔬 Enrichment Bar Plot


🚀 Features
📂 Upload ZIPs of BED files to generate consensus peaks

🧬 Peak annotation with ChIPseeker & TxDb.Hsapiens.UCSC.hg38.knownGene

🧠 Functional enrichment via enrichR (GO, KEGG, Reactome)

📊 Motif enrichment using motifmatchr + JASPAR2020

🧪 Differential accessibility analysis (DAA) via DESeq2

📈 Exploratory plots: PCA, UMAP, and heatmaps (plotly)

🤖 Machine learning: random forest classifier (AUC, feature importance)

🔍 Power analysis for experimental design

💾 Download buttons for all result tables

🎨 Stylish glowy anime-inspired UI

🖥️ Compatible with Docker, HPC (Singularity), or local RStudio

🧪 Sample Data
A sample ATAC-seq BED file is included:

📁 sample_data/ENCFF002CUU.bed

📚 Source: ENCODE GM12878 (hg19)

🔗 ENCFF002CUU

Use it to explore features before uploading your own.

📦 Requirements (if running outside Docker)
 in rstudio's control.

install.packages(c("shiny", "shinyjs", "plotly", "DT", "enrichR", "randomForest", "pROC", "powerAnalysis"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "ChIPseeker", "GenomicRanges", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db", "motifmatchr", "BSgenome.Hsapiens.UCSC.hg38",
  "JASPAR2020", "DESeq2"
))
🧬 How to Run the App
🔁 Option 1: Docker (Recommended)

in bash
git clone https://github.com/your-user/ATAC_APP.git
cd ATAC_APP
bash run.sh
Open in browser:
http://localhost:8787

💻 Option 2: Local RStudio
 in rstudio console

setwd("path/to/ATAC_APP")
library(shiny)
runApp(".")
🧠 Option 3: HPC with Singularity
 in bash

singularity build atac-shiny.sif Singularity.def
singularity run --bind $(pwd):/mnt atac-shiny.sif
Then forward ports and open in browser:
using 
ssh -L 8080:localhost:8080 user@cluster
http://localhost:8080
🎛️ How to Use the App
🔹 Step 1: Consensus Peaks
Upload a ZIP of BED files

Click Make Consensus Peaks

🔹 Step 2: Peak Annotation
Click Run Peak Annotation

View results in:

Peak Annotation Table

Annotation Pie Chart

🔹 Step 3: Enrichment Analysis
Select a database

Click Run Enrichment Analysis

View results in:

Enrichment Table

Enrichment Bar Plot

🔹 Step 4: Motif Enrichment
Select JASPAR family

Click Run Motif Enrichment

View motif match table and bar plot

🔹 Step 5: Counts & Metadata
Upload count matrix + metadata

Click Run DAA

Visualize differential peaks

🔹 Step 6: Dimensionality Reduction
PCA → Run PCA

UMAP → Run UMAP

Heatmap → Plot Heatmap

🔹 Step 7: Machine Learning
Click Run Random Forest

View:

RF Metrics

ROC Curve

Feature Importance Plot

🔹 Step 8: Power Analysis
Set effect size, FDR, replicates

Click Run Power Analysis

View:

Power Plot (w/ CI)

Power Table

📊 Tab Overview
Tab	Description
Peak Annotation Table	Annotated peaks with download
Annotation Pie Chart	Distribution of peak locations
Enrichment Table	enrichR results table
Enrichment Bar Plot	Top 10 enriched terms plot
Motif Enrichment Table	TF motif matches
Motif Enrichment Plot	Barplot of top motifs
Uploaded DAA Results	DESeq2 output from uploaded counts
PCA Plot	Principal component analysis
UMAP Plot	UMAP projection
Heatmap	Log-scaled count heatmap
RF Metrics	AUC, accuracy, and confusion matrix
Feature Importance	Top peaks by Gini importance
ROC Curve	ROC performance plot
Power Plot	Estimated power vs. replicates
Power Table	Power estimates with confidence intervals
README	This documentation inside the app

🛠️ Developer Notes
🔒 Error Logging
Errors logged to error_log.txt

email_log.R can be cron-scheduled to:

Email error summaries

Archive old logs to /logs/

📁 Folder Structure
bash
Copy
Edit
ATAC_APP/
├── app.R                  # Main Shiny app
├── Dockerfile             # Container build file
├── run.sh                 # Launch script
├── email_log.R            # Log emailer
├── error_log.txt          # Error logs
├── www/
│   └── fairy_tail.css     # Themed UI
├── logs/                  # Archived logs
├── sample_data/
│   └── ENCFF002CUU.bed    # Example ATAC peaks
└── README.md              # This file
👨‍🔬 Citation / Credit
If you use this tool in a publication, please consider citing the repo or giving a GitHub star ⭐.
MIT License — feel free to fork, remix, and enhance.

🧠 FAQ
Q: Can I use this with BAM files?
No, this app operates from BED/peak-level data. For BAM-to-peak workflows, use command-line tools (e.g., MACS2, deepTools) first.

Q: Why not cloud-deployed?
Due to reliance on Bioconductor packages and genome indices, deployment to serverless cloud platforms is non-trivial. We recommend Docker or HPC.

Q: How big can my datasets be?
Typical use cases scale well with up to ~10k peaks × 100 samples. Beyond that, RAM may become a constraint.

🧬 Built by scientists, for scientists — to make ATAC-seq suck less.
