FROM rocker/r-ver:4.0.4

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="build-essential apt-utils libcurl4-openssl-dev libssl-dev libxml2-dev libhdf5-dev libglpk-dev" && \
    BUILDPKGS_python="libncurses-dev libgdbm-dev libgdbm-compat-dev libz-dev tk-dev libsqlite3-dev libreadline-dev liblzma-dev libffi-dev libssl-dev libbz2-dev curl" && \
    apt-get update && \
    apt-get install -y $BUILDPKGS && \
    apt-get install -y $BUILDPKGS_python

# Install python from source
ENV PYVER 3.9.4
ENV PYVER_MAJOR 3.9
RUN mkdir /py_build && cd /py_build && \
    curl -O https://www.python.org/ftp/python/${PYVER}/Python-${PYVER}.tar.xz && \
    tar -xf Python-${PYVER}.tar.xz && cd Python-${PYVER} && \
    ./configure --enable-optimizations --enable-shared --prefix=/opt/python${PYVER} LDFLAGS=-Wl,-rpath=/opt/python${PYVER}/lib && \
    make -j 20 && \
    make altinstall && \
    ln -s /opt/python${PYVER}/bin/python${PYVER_MAJOR} /usr/local/bin/python && \
    ln -s /opt/python${PYVER}/bin/pip${PYVER_MAJOR} /usr/local/bin/pip && \
    ln -s /opt/python${PYVER}/bin/easy_install-${PYVER_MAJOR} /usr/local/bin/easy_install && \
    rm -rf /py_build
    
# Add package dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pandas && \
    pip install --no-cache-dir leidenalg && \
    pip install --no-cache-dir umap-learn && \
    R -e "install.packages('optparse')" && \
    R -e "install.packages('devtools')" && \
    R -e "install.packages('remotes')" && \
    R -e "install.packages('openxlsx')" && \
    R -e "install.packages('rmarkdown')" && \
    R -e "install.packages('BiocManager')" && \
    R -e "remotes::install_version('Seurat', '4.0.1')" && \
    R -e "devtools::install_github(repo = 'aertslab/SCopeLoomR')" && \
    R -e "BiocManager::install('DESeq2')" && \
    R -e "BiocManager::install('MAST')"

# Install last required packages + cleanup
RUN apt-get -y --no-install-recommends install \
        procps \
        pandoc && \
    apt-get remove --purge -y $BUILDPKGS && \
    apt-get remove --purge -y $BUILDPKGS_python && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean

# TODO: Flt-SNE (https://github.com/KlugerLab/FIt-SNE)