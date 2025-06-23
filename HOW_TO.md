🧬 JCAP_ATAC_SEQ APP — Usage Guide
🧬 Step 1: Upload Your Data
Upload ZIP of BED files for consensus peak calling.

Click "Make Consensus Peaks".

✅ You’ll receive a confirmation with the number of merged peaks.

📍 Step 2: Peak Annotation
Click "Run Peak Annotation".

Uses ChIPseeker with TxDb.Hsapiens.UCSC.hg38.knownGene.

Annotates peaks with genomic features.

Inspect results in:

📊 Peak Annotation Table — searchable & downloadable.

🥧 Annotation Pie Chart — visual distribution.

🧠 Step 3: Enrichment Analysis
Choose a database: GO, KEGG, or Reactome from the dropdown.

Click "Run Enrichment Analysis".

Uses enrichR based on annotated genes.

Inspect results in:

📋 Enrichment Table — full term list.

📈 Enrichment Bar Plot — top 10 terms by adjusted p-value.

🔍 Step 4: Motif Enrichment (Optional)
Select JASPAR family (e.g. ALL, TP53) from dropdown.

Click "Run Motif Enrichment".

Uses motifmatchr, JASPAR2020, and BSgenome.Hsapiens.UCSC.hg38.

Inspect results in:

📊 Motif Enrichment Table

📈 Top Motif Plot

🧪 Step 5: Counts & Metadata
Upload a CSV count matrix (peaks × samples).

Upload a CSV metadata file (sample × condition).

Click "Run DAA on Uploaded Counts".

Uses DESeq2 for differential analysis.

Inspect results in:

📉 Uploaded DAA Results

📈 Step 6: Exploratory Visualizations
Click "Run PCA" → PCA Plot

Click "Run UMAP" → UMAP Plot

Click "Plot Heatmap" → Heatmap (log2(count+1))

🤖 Step 7: Machine Learning — Random Forest
Click "Run Random Forest Classifier".

Trains on uploaded counts with 2-class condition labels.

Inspect results in:

📋 RF Metrics (AUC, accuracy)

📈 ROC Curve

🔥 Feature Importance Plot

🔬 Step 8: Power Analysis
Set Effect Size (log2FC), FDR α, and Replicates per Group.

Click "Run Power Analysis".

Inspect results in:

📈 Power Plot (with 95% confidence ribbons)

📊 Power Table (replicates, power, CI)

💾 Downloadable CSV via the “Download Power Table” button

💾 Downloadable Results
Use download buttons to save:

Peak Annotation Table

Enrichment Table

Motif Table

DAA Results

Power Table

📂 Tab Summary
Tab Name	Description
README	Usage instructions (this doc)
Peak Annotation Table	Annotated peaks with download
Annotation Pie Chart	Genomic feature distribution
Motif Enrichment Table	Matched TF motifs
Motif Enrichment Plot	Top motif bar chart
Enrichment Table	enrichR results
Enrichment Bar Plot	Top enriched terms plot
Uploaded DAA Results	DESeq2 differential analysis
PCA Plot	Sample PCA visualization
UMAP Plot	UMAP visualization of samples
Heatmap	Log2(transformed) heatmap of counts
RF Metrics	AUC & accuracy of classifier
Feature Importance	Top important peaks in RF
ROC Curve	ROC performance curve
Power Plot	Power vs replicates visualization
Power Table	Numeric power estimates with CI

⚠️ Errors?
All errors are logged to error_log.txt.
Admins may review and manage logs via email_log.R.

💡 Tips
Feel free to rerun any step without reuploading inputs.

Start with example datasets (if unfamiliar) to ensure proper format.

Hit “Run Power Analysis” to see how well your current experiment would detect effects — or plan future experiments accordingly.

Built by scientists, for scientists 🧬 — robust, reproducible ATAC-seq analysis in Shiny, no coding required.
