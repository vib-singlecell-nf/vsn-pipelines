FROM vibsinglecellnf/samtools:base-0.3

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="git \
        autoconf \
        automake \
        perl \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        bedtools \
        libncurses5-dev" && \
    apt-get update && \
    apt-get upgrade -y --no-install-recommends && \
    apt-get install -y --no-install-recommends $BUILDPKGS

# install htslib
ENV HTSLIB_VERSION 1.15.1
RUN curl -L -o /tmp/htslib-${HTSLIB_VERSION}.tar.bz2 \
        https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    mkdir -p /tmp/htslib-${HTSLIB_VERSION} && \
    tar jxvf /tmp/htslib-${HTSLIB_VERSION}.tar.bz2 -C /tmp/htslib-${HTSLIB_VERSION} --strip-components 1 && \
    cd /tmp/htslib-${HTSLIB_VERSION} && \
    ./configure \
        CFLAGS="-fPIC" && \
    make && \
    make install && \
    cd .. && rm -rf htslib-${HTSLIB_VERSION}*

# install samtools
ENV SAMTOOLS_VERSION 1.15.1
RUN curl -L -o /tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
        https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    mkdir -p /tmp/samtools-${SAMTOOLS_VERSION} && \
    tar jxvf /tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /tmp/samtools-${SAMTOOLS_VERSION} --strip-components 1 && \
    cd /tmp/samtools-${SAMTOOLS_VERSION} && \
    ./configure \
        --with-htslib=system && \
    make && \
    make install && \
    cd .. && rm -rf samtools-${SAMTOOLS_VERSION}*

RUN rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/* && \
    ldconfig

