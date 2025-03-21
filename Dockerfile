FROM rocker/r-ver:4.3.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
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
    libpng-dev \
    libjpeg-dev \
    libcairo2-dev \
    curl \
    pandoc \
    && apt-get clean

# Install R packages
RUN R -e "install.packages(c('shiny', 'shinyjs', 'plotly', 'DT', 'enrichR'))"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('ChIPseeker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db', 'GenomicRanges', 'clusterProfiler'), dependencies=TRUE)"

# Copy repo contents
COPY . /app
WORKDIR /app

# Port exposure
EXPOSE 8080


CMD ["R", "-e", "cat('🚀 Starting on port ', Sys.getenv('PORT'), '\\n'); port <- as.integer(Sys.getenv('PORT', unset = '8080')); options(shiny.port=port, shiny.host='0.0.0.0', shiny.launch.browser=FALSE); shiny::runApp('/app')"]
