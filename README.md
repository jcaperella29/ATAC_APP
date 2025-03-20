# 🔬 ATAC-seq Peak Annotation & Enrichment Viewer

This interactive Shiny app allows you to upload MACS2 `.narrowPeak` files, annotate peaks using **ChIPseeker**, and visualize functional enrichment with **enrichR**.

💡 Built with love in R, glowy CSS, and battle-tested with real-world ATAC-seq data.

---

## 🚀 Features

- 📂 Upload MACS2 `.narrowPeak` files
- 🧬 Annotate genomic regions (via `ChIPseeker`)
- 📊 Visualize annotation results (pie chart, table)
- 🧠 Functional enrichment with GO, KEGG, Reactome (via `enrichR`)
- ⚠️ Designed with stability and reliability in mind

---

## 📦 Requirements

Install the required R packages:

```r
install.packages(c(
  "shiny", "shinyjs", "plotly", "DT", "enrichR", "clusterProfiler",
  "ChIPseeker", "GenomicRanges", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db", "blastula"
))


🧪 How to Run
Clone the repo and launch the app locally:
in your terminal run:
R

library(shiny)
runApp("path/to/ATAC_APP"

## 🧪 How to Use the App

1. **Upload MACS2 narrowPeak file**  
   - Click **"Upload MACS2 narrowPeak File"**  
   - Select your `.narrowPeak` file to load peaks into the app

2. **Run Peak Annotation**  
   - Click the **"Run Peak Annotation"** button  
   - Uses `ChIPseeker` to assign each peak to a genomic feature (e.g., promoter, intron, intergenic)

3. **Run Enrichment Analysis**  
   - Use the dropdown to choose a database (GO, KEGG, Reactome)  
   - Then click **"Run Enrichment Analysis"**  
   - Identifies overrepresented biological processes or pathways using your annotated gene list

---

## 📊 Tabs Overview

- **Peak Annotation Table**  
  View all annotated peaks in a searchable, filterable table  
  → Includes download button for `.csv` export

- **Annotation Pie Chart**  
  Visual summary of the distribution of genomic features assigned to peaks

- **Enrichment Table**  
  Tabular view of enriched terms from the selected database  
  → Also downloadable as `.csv`

- **Enrichment Bar Plot**  
  Visualize the top 10 enriched terms by adjusted p-value (interactive `plotly` bar chart)

- **README**  
  You're here! App walkthrough and feature guide embedded directly for convenience.


🛠️ Developer Notes
🔒 Internal Error Logging & Monitoring
The app logs all server-side errors to error_log.txt.
If deployed to a server, the script email_log.R can be scheduled (e.g. with cron or taskscheduleR) to:

Email errors to the developer automatically at midnight
Archive old logs into /logs/
This system is invisible to users and intended for maintainers only.

📁 Folder Structure
bash
Copy
Edit
ATAC_APP/
├── app.R
├── email_log.R            # Internal only (dev monitoring)
├── www/
│   └── fairy_tail.css     # Glowy anime-themed UI
├── error_log.txt          # Generated automatically
├── logs/                  # Archived logs (auto-generated)
├── .gitignore
└── README.md
