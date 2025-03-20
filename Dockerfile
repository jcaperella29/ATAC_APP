FROM rocker/shiny:4.3.1

# Install system libs
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
    pandoc \
    libpng-dev \
    libjpeg-dev \
    libcairo2-dev \
    && apt-get clean

# Install R packages
RUN R -e "install.packages(c('shiny', 'shinyjs', 'plotly', 'DT', 'enrichR'))"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('ChIPseeker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db', 'GenomicRanges', 'clusterProfiler'), dependencies=TRUE)"

# Copy app
COPY . /srv/shiny-server/

# Set permissions
RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 8080

# Cloud Run health check (optional but good)
HEALTHCHECK CMD curl --fail http://localhost:8080 || exit 1

# CMD is already set to shiny-server in base image
