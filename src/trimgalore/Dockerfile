FROM python:3.7-slim

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="build-essential zlib1g-dev git curl" && \
    apt-get update && \
    apt-get install -y --no-install-recommends apt-utils debconf locales && dpkg-reconfigure locales && \
    apt-get install -y --no-install-recommends $BUILDPKGS

RUN pip install -U pip

##################################################
# cutadapt
RUN pip install cutadapt

##################################################
# fastQC
# RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip

##################################################
# trim galore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz && \
    mv TrimGalore-0.6.6/trim_galore /usr/bin/ && \
    rm -r TrimGalore-0.6.6


RUN apt-get -y update && \
    apt-get -y --no-install-recommends install \
        # Need to run ps
        procps \
        pigz \
        less && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

