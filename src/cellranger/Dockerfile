FROM debian:buster-slim

ENV CELLRANGER_VER 5.0.1

# Pre-downloaded files:
# cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
COPY cellranger-$CELLRANGER_VER.tar.gz /opt/cellranger-$CELLRANGER_VER.tar.gz
# bcl2fastq2: https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html
COPY bcl2fastq2-v2-20-0-linux-x86-64.zip /tmp/bcl2fastq2-v2-20-0-linux-x86-64.zip

RUN apt-get update && \
    # Need to run ps
    apt-get -y install procps && \
    apt-get install -y bsdtar p7zip-full cpio wget unzip

# Install bcl2fastq by extracting rpm
RUN  \
    cd /tmp && \
    unzip bcl2fastq2-v2-20-0-linux-x86-64.zip && \
    7z e bcl2fastq2-v2.20.0*-Linux-x86_64.rpm && \
    cpio -idmv "./usr/local/bin/bcl2fastq" < bcl2fastq2-v2.20.0.422-1.x86_64.cpio && \
    mv usr/local/bin/bcl2fastq /usr/bin && \
    rm bcl2fastq2*

# Install CellRanger from tgz file
RUN \
  cd /opt && \
  bsdtar -xzvf cellranger-$CELLRANGER_VER.tar.gz && \
  export PATH=/opt/cellranger-$CELLRANGER_VER:$PATH && \
  ln -s /opt/cellranger-$CELLRANGER_VER/cellranger /usr/bin/cellranger && \
  rm -rf /opt/cellranger-$CELLRANGER_VER.tar.gz
