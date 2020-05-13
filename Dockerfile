FROM dweemx/sctx-seurat:3.1.2

RUN apt-get -y update && \
    apt-get install -y libcurl4-openssl-dev libxml2-dev zlib1g-dev libhdf5-dev && \
    apt-get install -y libssl-dev && \
    # png.h: No such file or directory
    apt-get install -y libpng-dev && \
    R -e "install.packages('doFuture')" && \
    R -e "install.packages('doRNG')" && \
    R -e "install.packages('optparse')" && \
    R -e "install.packages('dismo')" && \
    R -e "devtools::install_github(repo = 'aertslab/SCopeLoomR')" && \
    # Need to run ps
    apt-get -y install procps && \
    apt-get -y install libxml2 && \
    # Clean
    rm -rf /tmp/* && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean
