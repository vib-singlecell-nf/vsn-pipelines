FROM vibsinglecellnf/samtools:base-0.3

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="perl python3 python3-pip" && \
    apt-get update && \
    apt-get upgrade -y --no-install-recommends && \
    apt-get install -y --no-install-recommends $BUILDPKGS && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3.9 100 && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

# Update pip.
RUN pip install --no-cache-dir --upgrade pip

# Install cutadapt.
RUN pip install --no-cache-dir cutadapt==4.1

# Install TrimGalore.
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.7.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz && \
    mv TrimGalore-0.6.7/trim_galore /usr/bin/ && \
    rm -r TrimGalore-0.6.7

