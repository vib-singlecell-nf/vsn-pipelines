FROM python:3.8-slim

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="build-essential \
        samtools \
        tabix \
        git" && \
    apt-get update && \
    apt-get install -y --no-install-recommends apt-utils debconf locales && dpkg-reconfigure locales && \
    apt-get install -y --no-install-recommends $BUILDPKGS

RUN pip install --no-cache-dir -U pip

RUN git clone https://github.com/timoast/sinto.git /tmp/sinto && \
    cd /tmp/sinto && \
    pip install . && \
    cd .. && rm -r sinto

RUN apt-get -y update && \
    apt-get -y --no-install-recommends install \
        # Need to run ps
        procps \
        less && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

