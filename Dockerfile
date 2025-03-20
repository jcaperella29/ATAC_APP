FROM rocker/shiny:4.3.1

# Avoid prompts
ENV DEBIAN_FRONTEND=noninteractive

# System dependencies for your R packages
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
    pandoc \
    && apt-get clean

# Install BiocManager and required R packages
RUN R -e "install.packages(c('shiny', 'shinyjs', 'plotly', 'DT', 'enrichR'))"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('ChIPseeker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db', 'GenomicRanges', 'clusterProfiler'))"

# Copy all files into image
COPY . /srv/shiny-server/

# Permissions
RUN chown -R shiny:shiny /srv/shiny-server

# Expose default port
EXPOSE 3838

# Run as `shiny` user
USER shiny

# Launch app
CMD ["/usr/bin/shiny-server"]
