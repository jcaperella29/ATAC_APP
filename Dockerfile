FROM rocker/r-ver:4.3.1

ENV DEBIAN_FRONTEND=noninteractive
ENV MAKEFLAGS="-j4"

# System packages
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libglu1-mesa \
    libxkbcommon-x11-0 \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libbz2-dev \
    liblzma-dev \
    libz-dev \
    libncurses-dev \
    libxt-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libicu-dev \
    zlib1g-dev \
    pandoc && apt-get clean

# Install core R packages
RUN R -e "install.packages(c('shiny', 'shinyjs', 'plotly', 'DT', 'enrichR'), repos='https://cloud.r-project.org')"

# Bioconductor with explicit enrichplot to avoid lazy load error
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(version='3.18', ask=FALSE)"
RUN R -e "BiocManager::install(c('ChIPseeker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db', 'GenomicRanges', 'clusterProfiler', 'enrichplot'), ask=FALSE, update=FALSE, dependencies=TRUE)"

# Copy app into image
COPY . /app

# Set working directory
WORKDIR /app

# Expose port expected by Cloud Run
EXPOSE 8080

# Run the app manually on port 8080
CMD ["Rscript", "-e", "cat('🔥 Booting on port=', Sys.getenv('PORT'), '\\n'); options(shiny.port=as.integer(Sys.getenv('PORT')), shiny.host='0.0.0.0'); shiny::runApp('/app')"]
