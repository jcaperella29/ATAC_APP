FROM rocker/r-ver:4.3.1
CMD ["R", "-e", "cat('🔥 App booting...\n'); options(shiny.port=as.integer(Sys.getenv('PORT')), shiny.host='0.0.0.0'); shiny::runApp('/app')"]

# System deps
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    libglu1-mesa libxkbcommon-x11-0 libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev pandoc && apt-get clean

# Install R packages
RUN R -e "install.packages(c('shiny', 'shinyjs', 'plotly', 'DT', 'enrichR'))"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('ChIPseeker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db', 'GenomicRanges', 'clusterProfiler'))"

# Copy app
COPY . /app

# Expose port
EXPOSE 8080

# Launch app manually using shiny::runApp

CMD ["R", "-e", "cat('🚀 Starting App on PORT=', Sys.getenv('PORT'), '\\n'); options(shiny.port=as.integer(Sys.getenv('PORT')), shiny.host='0.0.0.0'); shiny::runApp('/app')"]
