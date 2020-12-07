FROM rocker/r-ver:3.6.0

RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev libxml2-dev zlib1g-dev libhdf5-dev && \
    R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) { install.packages('BiocManager') }; BiocManager::install('DropletUtils')" \
    R -e "install.packages('dplyr')" \
    R -e "install.packages('ggplot2')" \
    # reticulate seems to be needed
    R -e "install.packages('reticulate')" \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*
