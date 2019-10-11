FROM python:3.7.4-slim-stretch AS compile-image

RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential gcc apt-utils cmake openssh-client git && \
    apt-get install -y python-dev libhdf5-dev libxml2-dev zlib1g-dev # Needed for igraph && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

RUN python -m venv /opt/venv
# Make sure we use the virtualenv:
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --upgrade pip && \
    pip install --no-cache-dir --upgrade numpy && \
    pip install --no-cache-dir seaborn scikit-learn statsmodels tables && \
    pip install --no-cache-dir python-igraph louvain && \
    pip install --no-cache-dir MulticoreTSNE && \
    pip install --no-cache-dir numba==0.43.1 && \
    pip install --no-cache-dir scanpy && \
    pip install --no-cache-dir mnnpy && \
    pip install --no-cache-dir annoy==1.15.2 && \
    pip install --no-cache-dir bbknn && \
    pip install --no-cache-dir loompy && \
    python3 -m pip install ipykernel && \
    pip install --no-cache-dir papermill

FROM python:3.7.4-slim-stretch AS build-image
RUN apt-get -y update && \
    # Need to run ps
    apt-get -y install procps && \
    apt-get -y install libxml2 && \
    # Need to run MulticoreTSNE
    apt-get -y install libgomp1 && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

COPY --from=compile-image /opt/venv /opt/venv

# Make sure we use the virtualenv:
ENV PATH="/opt/venv/bin:$PATH"