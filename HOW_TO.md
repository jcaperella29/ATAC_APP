# 🧬 JCAP_ATAC_SEQ APP — Usage Guide

## Step 1: Upload Your Data
- **Upload** a ZIP of BED files for consensus peak calling.
- Click **"Make Consensus Peaks"**.
- ✅ Confirmation message displays the number of merged peaks.

## Step 2: Peak Annotation
- Click **"Run Peak Annotation"**.
- Peaks are annotated with nearest gene and distance to TSS using `TxDb.Hsapiens.UCSC.hg38.knownGene` and `org.Hs.eg.db`.
- View results in:
    - **📊 Peak Annotation Table** — searchable & downloadable
    - **🥧 Annotation Pie Chart** — visual gene distribution

## Step 3: Motif Enrichment (Optional)
- Select JASPAR TF family (e.g. ALL, TP53) from dropdown.
- Click **"Run Motif Enrichment"**.
- Uses `motifmatchr`, `JASPAR2020`, and `BSgenome.Hsapiens.UCSC.hg38` for motif analysis.
- View results in:
    - **📊 Motif Enrichment Table**
    - **📈 Top Motif Plot**

## Step 4: Counts & Metadata
- **Upload** a CSV count matrix (peaks × samples).
- **Upload** a CSV metadata file (sample × condition).
- Click **"Run DAA on Uploaded Counts"**.
- Runs **DESeq2** for differential analysis.
- View results in: **📉 Uploaded DAA Results**

## Step 5: Exploratory Visualizations
- Click **"Run PCA"** → View PCA Plot
- Click **"Run UMAP"** → View UMAP Plot
- Click **"Plot Heatmap"** → Heatmap of log2(count+1)

## Step 6: Machine Learning — Random Forest
- Click **"Run Random Forest Classifier"**.
- Trains classifier on uploaded counts with 2-class condition labels.
- View results in:
    - **📋 RF Metrics (AUC, accuracy)**
    - **📈 ROC Curve**
    - **🔥 Feature Importance Plot**

## Step 7: Power Analysis
- Set **Effect Size** (log2FC), **FDR α**, and **Replicates per Group**.
- Click **"Run Power Analysis"**.
- View results in:
    - **📈 Power Plot** (with 95% confidence intervals)
    - **📊 Power Table** (replicates, power, CI)
    - Downloadable CSV via the “Download Power Table” button

## 💾 Downloadable Results
Use download buttons to save:
- Peak Annotation Table
- Motif Table
- DAA Results
- Power Table

---

## 📂 Tab Summary

| Tab Name                 | Description                             |
|--------------------------|-----------------------------------------|
| README                   | Usage instructions (this doc)           |
| Peak Annotation Table    | Annotated peaks with download           |
| Annotation Pie Chart     | Gene distribution pie chart             |
| Motif Enrichment Table   | Matched TF motifs (motifmatchr)         |
| Motif Enrichment Plot    | Top motif bar chart                     |
| Uploaded DAA Results     | DESeq2 differential analysis            |
| PCA Plot                 | Sample PCA visualization                |
| UMAP Plot                | UMAP visualization of samples           |
| Heatmap                  | Log2(count+1) heatmap                   |
| RF Metrics               | AUC & accuracy of classifier            |
| Feature Importance       | Top important peaks in RF               |
| ROC Curve                | ROC performance curve                   |
| Power Plot               | Power vs replicates visualization       |
| Power Table              | Numeric power estimates with CI         |

---

## ⚠️ Errors?
- All errors are logged to `error_log.txt`.
- Admins can review and manage logs via `email_log.R`.

## 💡 Tips
- You can rerun any step without reuploading inputs.
- Use included example datasets to verify input formats.
- “Run Power Analysis” helps you evaluate or plan experiments for statistical power.
- Built by scientists, for scientists 🧬 — robust, reproducible ATAC-seq analysis in Shiny, no coding required.

---

*Annotation and motif analysis are performed entirely with open and actively maintained Bioconductor and CRAN packages. No dependency on ChIPseeker or enrichR.*

