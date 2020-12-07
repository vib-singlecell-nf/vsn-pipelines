FROM debian:buster-slim

ENV CELLRANGERATAC_VER 1.2.0

# Pre-downloaded files:
# cellranger atac: https://support.10xgenomics.com/single-cell-atac/software/downloads/latest
COPY cellranger-atac-$CELLRANGERATAC_VER.tar.gz /opt/cellranger-atac-$CELLRANGERATAC_VER.tar.gz
# bcl2fastq2: https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html
COPY bcl2fastq2-v2-20-0-linux-x86-64.zip /tmp/bcl2fastq2-v2-20-0-linux-x86-64.zip

RUN apt-get update && \
    apt-get install -y bsdtar p7zip-full cpio wget unzip procps

# Install bcl2fastq by extracting rpm
RUN  \
    cd /tmp && \
    unzip bcl2fastq2-v2-20-0-linux-x86-64.zip && \
    7z e bcl2fastq2-v2.20.0*-Linux-x86_64.rpm && \
    cpio -idmv "./usr/local/bin/bcl2fastq" < bcl2fastq2-v2.20.0.422-1.x86_64.cpio && \
    mv usr/local/bin/bcl2fastq /usr/bin && \
    rm bcl2fastq2*

# Install Cell Ranger ATAC from tgz file
RUN \
  cd /opt && \
  bsdtar -xzvf cellranger-atac-$CELLRANGERATAC_VER.tar.gz && \
  export PATH=/opt/cellranger-atac-$CELLRANGERATAC_VER:$PATH && \
  ln -s /opt/cellranger-atac-$CELLRANGERATAC_VER/cellranger-atac /usr/bin/cellranger-atac && \
  rm -rf /opt/cellranger-atac-$CELLRANGERATAC_VER.tar.gz

