FROM python:3.7.4-slim AS compile-image

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="build-essential apt-utils \
    python3-dev libhdf5-dev libfreetype6-dev libtool \
    m4 autoconf automake patch bison flex libpng-dev libopenblas-dev \
    tcl-dev tk-dev libxml2-dev zlib1g-dev libffi-dev cmake" && \
    apt-get update && \
    apt-get install -y --no-install-recommends apt-utils debconf locales && dpkg-reconfigure locales && \
    apt-get install -y --no-install-recommends $BUILDPKGS

RUN python -m venv /opt/venv
# Make sure we use the virtualenv:
ENV PATH="/opt/venv/bin:$PATH"

# install dependencies:
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /tmp/requirements.txt

RUN pip install --no-cache-dir scanpy==1.8.1


FROM python:3.7.4-slim AS build-image

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && \
    apt-get -y --no-install-recommends install \
    # Need to run ps
    procps \
    libxml2 \
    # Need to run MulticoreTSNE
    libgomp1 && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

COPY --from=compile-image /opt/venv /opt/venv
# Make sure we use the virtualenv:
ENV PATH="/opt/venv/bin:$PATH"

