FROM vibsinglecellnf/samtools:1.11
# source: https://github.com/vib-singlecell-nf/bwamaptools/blob/master/Dockerfile.samtools

RUN apk add --virtual build-dependencies git

# install bwa
RUN git clone --depth 1 https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    mv /bwa/bwa /usr/local/bin/

RUN apk del build-dependencies && \
    rm -rf /var/cache/apk/*

