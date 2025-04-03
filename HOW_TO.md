🧬 Step 1: Upload Your Data
Click "Upload MACS2 narrowPeak File or BED file "
Select a .narrowPeak file from your MACS2 output or a .bed file
Wait for confirmation message
✅ File Accepted → You're ready to annotate.

📍 Step 2: Run Peak Annotation
Click "Run Peak Annotation"
Internally uses ChIPseeker + TxDb.Hsapiens.UCSC.hg38.knownGene
Annotates peaks with genomic features (e.g., promoter, exon, intergenic)
Check:

📊 Peak Annotation Table tab → full searchable table
🥧 Annotation Pie Chart tab → visual breakdown
🧠 Step 3: Run Enrichment Analysis
Select database (GO, KEGG, or Reactome)
Click "Run Enrichment Analysis"
Uses enrichR to analyze genes from the annotation step
Check:

📋 Enrichment Table → full term results
📈 Enrichment Bar Plot → top 10 enriched terms by adjusted p-value
💾 Step 4: Download Results
Look for Download CSV buttons in:

Peak Annotation Table
Enrichment Table
Top Enrichment Bar Plot (data)
📂 Tabs Overview
Tab	Description
Peak Annotation Table	View and download annotated peaks
Annotation Pie Chart	Pie chart of genomic region distribution (e.g., promoters, exons, etc.)
Enrichment Table	Enrichr-based table of enriched biological terms
Enrichment Bar Plot	Plotly bar chart of top enriched terms
README	You're here — usage guide embedded inside the app
⚠️ Errors?
Any issues are logged automatically to error_log.txt.
For admins: use email_log.R to send daily error emails.

💡 Tip:
You can rerun enrichment with different databases anytime — no need to re-upload or re-annotate unless your data changes.

Built by a scientist, for scientists 🧬
Fast, clean, reproducible ATAC-seq annotation — with zero coding required.
