FROM vibsinglecellnf/samtools:0.3-1.15.1

ENV DEBIAN_FRONTEND=noninteractive

# install bwa
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    mv /bwa/bwa /usr/local/bin/

# install bwa-mem2
ENV BWAMEM2_VER 2.2.1
RUN cd /tmp && \
    curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWAMEM2_VER}/bwa-mem2-${BWAMEM2_VER}_x64-linux.tar.bz2 \
      | tar jxf - --no-same-owner && \
    mv bwa-mem2-${BWAMEM2_VER}_x64-linux/bwa-mem2* /usr/local/bin

RUN rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/* && \
    ldconfig

