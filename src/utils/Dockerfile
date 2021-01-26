FROM python:3.7.9-slim-stretch AS compile-image

RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential gcc apt-utils cmake openssh-client git && \
    apt-get install -y python-dev libxml2-dev zlib1g-dev && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

RUN python -m venv /opt/venv
# Make sure we use the virtualenv:
ENV PATH="/opt/venv/bin:$PATH"

RUN python3 -m pip install ipykernel && \
    pip install --no-cache-dir papermill && \
    pip install --no-cache-dir pysradb==1.0.0 && \
    pip install --no-cache-dir nbconvert==5.6.0 && \
    pip install --no-cache-dir nbmerge==0.0.4 && \
    pip install --no-cache-dir nbformat==5.0.8

FROM python:3.7.9-slim-stretch AS build-image
RUN apt-get -y update && \
    # Need to run ps
    apt-get -y install procps && \
    apt-get -y install libxml2 && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

COPY --from=compile-image /opt/venv /opt/venv

# Make sure we use the virtualenv:
ENV PATH="/opt/venv/bin:$PATH"