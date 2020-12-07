*** System restart required ***
Last login: Wed Nov 25 10:25:17 2020 from 10.118.229.21
    u0125489 @ gbw-s-seq10 : ~
$ tmux a
FROM ubuntu:20.04


ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev \
    git \
    libhts-dev \
    curl \
    wget" && \
    apt-get update && \
    apt-get install -y apt-utils debconf locales ca-certificates && dpkg-reconfigure locales && \
    apt-get install -y --no-install-recommends $BUILDPKGS && \
    apt-get install -y --no-install-recommends clang python3 llvm-6.0 pigz

# install htslib
ENV HTSLIB_VERSION 1.11
RUN curl -L -o /tmp/htslib-${HTSLIB_VERSION}.tar.bz2 \
        https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    mkdir -p /tmp/htslib-${HTSLIB_VERSION} && \
    tar jxvf /tmp/htslib-${HTSLIB_VERSION}.tar.bz2 -C /tmp/htslib-${HTSLIB_VERSION} --strip-components 1 && \
    cd /tmp/htslib-${HTSLIB_VERSION} && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf htslib-${HTSLIB_VERSION}

# install samtools
ENV SAMTOOLS_VERSION 1.11
RUN curl -L -o /tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
        https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    mkdir -p /tmp/samtools-${SAMTOOLS_VERSION} && \
    tar jxvf /tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /tmp/samtools-${SAMTOOLS_VERSION} --strip-components 1 && \
    cd /tmp/samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf samtools-${SAMTOOLS_VERSION}

# install seq (https://github.com/seq-lang/seq/):
ENV SEQ_VERSION=0.9.11
RUN mkdir -p /opt/seq && \
    wget https://github.com/seq-lang/seq/releases/download/v${SEQ_VERSION}/seq-linux-x86_64.tar.gz && \
    tar xzf seq-linux-x86_64.tar.gz --strip-components 1 -C /opt/seq && \
    rm seq-linux-x86_64.tar.gz
ENV PATH="/opt/seq/bin:${PATH}"
# ENV OMP_NUM_THREADS=1
ENV SEQ_PYTHON=/usr/lib/x86_64-linux-gnu/libpython3.8.so.1

# install single_cell_toolkit
# https://github.com/aertslab/single_cell_toolkit
RUN git clone --depth=1 https://github.com/aertslab/single_cell_toolkit.git /opt/single_cell_toolkit
ENV seq_root_dir=/opt/seq
ENV PATH="/opt/single_cell_toolkit:${PATH}"

# final set of packages for usability
RUN apt-get -y update && \
    apt-get -y --no-install-recommends install \
        procps \
        bash-completion \
        less && \
    apt-get remove --purge -y $BUILDPKGS && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

