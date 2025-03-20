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

r
Copy
Edit
library(shiny)
runApp("path/to/ATAC_APP"

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
