FROM alpine:3.12.1

ENV SAMTOOLS_VERSION 1.11
ENV HTSLIB_VERSION 1.11

RUN apk update && \
    apk add --update autoconf automake make gcc musl-dev perl bash zlib-dev bzip2-dev xz-dev curl-dev libressl-dev ncurses-dev mawk && \
    apk add --virtual build-dependencies curl git

RUN curl -L -o /tmp/htslib-${HTSLIB_VERSION}.tar.bz2 \
        https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    mkdir -p /tmp/htslib-${HTSLIB_VERSION} && \
    tar jxvf /tmp/htslib-${HTSLIB_VERSION}.tar.bz2 -C /tmp/htslib-${HTSLIB_VERSION} --strip-components 1 && \
    cd /tmp/htslib-${HTSLIB_VERSION} && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf htslib-${HTSLIB_VERSION}

RUN curl -L -o /tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
        https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    mkdir -p /tmp/samtools-${SAMTOOLS_VERSION} && \
    tar jxvf /tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /tmp/samtools-${SAMTOOLS_VERSION} --strip-components 1 && \
    cd /tmp/samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf samtools-${SAMTOOLS_VERSION}

RUN git clone https://github.com/aertslab/single_cell_toolkit.git /tmp/single_cell_toolkit && \
    mv /tmp/single_cell_toolkit/*sh /usr/bin

RUN apk del build-dependencies && \
    rm -rf /var/chache/apk/*

