FROM debian:buster-slim

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="build-essential \
        autoconf cmake git \
        automake make gcc perl bedtools \
        libbz2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        zlib1g-dev \
        liblzma-dev \
        curl \
        libncurses5-dev"&& \
    apt-get update && \
    apt-get install -y --no-install-recommends apt-utils debconf locales && dpkg-reconfigure locales && \
    apt-get install -y --reinstall ca-certificates && \
    apt-get install -y --no-install-recommends $BUILDPKGS

# Install htslib
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git && \
    cd htslib && \
    autoheader && \
    autoconf && \
    ./configure --prefix=/usr/local/ && \
    make && \
    make install

# Install SAMtools
RUN git clone https://github.com/samtools/samtools.git && \
    cd samtools && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure --prefix=/usr/local/ && \
    make && \
    make install

# install bwa
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    mv /bwa/bwa /usr/local/bin/

# install bwa-mem2
ENV BWAMEM2_VER 2.1
RUN cd /tmp && \
    curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWAMEM2_VER}/bwa-mem2-${BWAMEM2_VER}_x64-linux.tar.bz2 | tar jxf - && \
    mv bwa-mem2-${BWAMEM2_VER}_x64-linux/bwa-mem2* /usr/local/bin

RUN apt-get -y update && \
    apt-get -y --no-install-recommends install \
        # Need to run ps
        procps \
        less && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

