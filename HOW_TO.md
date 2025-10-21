🧬 JCAP_ATAC_SEQ APP — Multi-Species Usage Guide
Step 1 | Upload Your Data

Upload a ZIP of BED files for consensus peak calling.

Click “Make Consensus Peaks.”

✅ Confirmation message will display the number of merged peaks.

Step 2 | Peak Annotation (Organism-Aware)

Click “Run Peak Annotation.”

Peaks are annotated with nearest gene and distance to TSS using species-specific Bioconductor packages:

Species	TxDb Package	OrgDb Package
Human	TxDb.Hsapiens.UCSC.hg38.knownGene	org.Hs.eg.db
Mouse	TxDb.Mmusculus.UCSC.mm10.knownGene	org.Mm.eg.db
Zebrafish	TxDb.Drerio.UCSC.danRer11.refGene	org.Dr.eg.db
Fly	TxDb.Dmelanogaster.UCSC.dm6.ensGene	org.Dm.eg.db

View results in:

📊 Peak Annotation Table — searchable & downloadable

🥧 Annotation Pie Chart — visual gene distribution

Step 3 | Motif Enrichment (Optional)

Select the species and TF family (e.g., ALL, TP53) from the dropdown.

Click “Run Motif Enrichment.”

Uses:

TFBSTools + JASPAR2020 (species filtered via taxID)

motifmatchr for genome-wide PWM matching

BSgenome.* for the chosen organism:

BSgenome.Hsapiens.UCSC.hg38

BSgenome.Mmusculus.UCSC.mm10

BSgenome.Drerio.UCSC.danRer11

BSgenome.Dmelanogaster.UCSC.dm6

View results in:

📊 Motif Enrichment Table

📈 Top Motif Plot

Step 4 | Counts & Metadata

Upload a CSV count matrix (peaks × samples).

Upload a CSV metadata file (samples × condition).

Click “Run DAA on Uploaded Counts.”

Runs DESeq2 for differential accessibility analysis.

View results in 📉 Uploaded DAA Results.

Step 5 | Exploratory Visualizations

Run PCA → PCA plot of samples

Run UMAP → UMAP embedding of samples

Plot Heatmap → Heatmap of log₂(count + 1)

Step 6 | Machine Learning — Random Forest

Click “Run Random Forest Classifier.”

Trains a 2-class classifier on uploaded counts.

View results in:

📋 RF Metrics (AUC, Accuracy)

📈 ROC Curve

🔥 Feature Importance Plot

Step 7 | Power Analysis

Set Effect Size, FDR α, and Replicates per Group.

Click “Run Power Analysis.”

View results in:

📈 Power Plot (with 95 % CI)

📊 Power Table (replicates vs power)

Downloadable CSV via Download Power Table

💾 Downloadable Results

You can save:

Peak Annotation Table

Motif Enrichment Table

DAA Results

Power Table

📂 Tab Summary
Tab Name	Description
README	Usage instructions (this document)
Peak Annotation Table	Annotated peaks with download
Annotation Pie Chart	Gene distribution pie chart
Motif Enrichment Table	Matched TF motifs across species
Motif Enrichment Plot	Top motif bar chart
Uploaded DAA Results	DESeq2 differential analysis
PCA Plot	Principal component analysis
UMAP Plot	UMAP sample embedding
Heatmap	log₂(count + 1) heatmap
RF Metrics	AUC and accuracy
Feature Importance	Top important peaks
ROC Curve	ROC performance curve
Power Plot	Power vs replicates
Power Table	Numeric power estimates
⚠️ Error Logging

All runtime errors are written to error_log.txt.

Admin scripts (e.g., email_log.R) can aggregate or forward logs.

💡 Tips

Rerun any module without reuploading data.

Use included example datasets to verify formats.

“Run Power Analysis” helps plan ATAC-seq replicates for desired power.

Built by scientists, for scientists 🧬 — reproducible, species-aware ATAC-seq analytics with no coding required.
