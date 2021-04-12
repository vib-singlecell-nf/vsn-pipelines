FROM rocker/r-base:4.0.4

RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev libssl-dev libhdf5-dev && \
    apt-get -y install procps && \
    apt-get -y install libxml2 libxml2-dev && \
    # R packages
    R -e "install.packages('optparse')" && \
    R -e "install.packages('devtools')" && \
    R -e "install.packages('remotes')" && \
    R -e "remotes::install_version('Seurat', '4.0.1')" && \
    R -e "devtools::install_github(repo = 'aertslab/SCopeLoomR')" && \
    # Cleanup
    rm -rf /tmp/* && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean

# TODO: leidenalg
# TODO: DESeq2
# TODO: MAST
# TODO: Flt-SNE
# TODO: umap learn (python)