FROM rocker/r-ver:3.6.3

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="build-essential apt-utils libssl-dev libcurl4-openssl-dev zlib1g-dev libssh2-1-dev libv8-dev libgsl-dev libxml2-dev libpng-dev libunwind-dev libbz2-dev liblzma-dev libicu-dev libpcre3-dev libjpeg-dev libcairo2-dev libxt-dev libx11-dev libgit2-dev git" && \
    apt-get update && \
    apt-get install -y $BUILDPKGS

##################################################
# build python for MACS2
RUN BUILDPKGS_python="build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev wget curl" && \
    apt update && \
    apt install -y $BUILDPKGS_python

# compile python from source:
ENV PYVER 3.8.3
ENV PYVER_MAJOR 3.8
RUN mkdir /py_build && cd /py_build && \
    curl -O https://www.python.org/ftp/python/${PYVER}/Python-${PYVER}.tar.xz && \
    tar -xf Python-${PYVER}.tar.xz && cd Python-${PYVER} && \
    ./configure --enable-optimizations && \
    make -j 20 && \
    make altinstall && \
    ln -s /usr/local/bin/python${PYVER_MAJOR} /usr/local/bin/python && \
    ln -s /usr/local/bin/pip${PYVER_MAJOR} /usr/local/bin/pip && \
    rm -rf /py_build

# install numpy and MACS2
RUN git clone https://github.com/taoliu/MACS.git /MACS && \
    pip install --trusted-host pypi.python.org --upgrade pip && pip install --trusted-host pypi.python.org -r /MACS/requirements.txt && \
    cd /MACS && python setup.py install && \
    rm -r /MACS

##################################################
# prepare R install
RUN R -e "install.packages( c('devtools','BiocManager', 'pheatmap', 'XML', 'optparse') )"

# install ArchR itself
ARG commit=master
RUN R -e "devtools::install_github('GreenleafLab/ArchR', ref='$commit', repos = BiocManager::repositories())" && \
    R -e "library(ArchR); ArchR::installExtraPackages()"

# additional support packages:
RUN R -e "BiocManager::install(c('BiocVersion', 'BSgenome.Hsapiens.UCSC.hg38', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Mmusculus.UCSC.mm9', 'BSgenome.Mmusculus.UCSC.mm10'), update=TRUE, ask=FALSE)"

# final set of packages for usability
RUN apt-get -y --no-install-recommends install \
        procps \
        bash-completion \
        curl \
        tk \
        less && \
    apt-get remove --purge -y $BUILDPKGS && \
    apt-get remove --purge -y $BUILDPKGS_python && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

CMD ["R"]

