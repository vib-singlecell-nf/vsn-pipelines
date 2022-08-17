FROM vibsinglecellnf/samtools:0.3-1.15.1

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 \
        python3-venv \
        nasm \
        libtool \
        wget && \
        update-alternatives --install /usr/bin/python python /usr/bin/python3.9 100

# install igzip (https://github.com/intel/isa-l)
RUN git clone --depth=1 https://github.com/intel/isa-l.git /tmp/isa-l && \
    cd /tmp/isa-l && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -r isa-l


RUN python -m venv /opt/venv
# Make sure we use the virtualenv:
ENV PATH="/opt/venv/bin:$PATH"

ENV POLARS_MAX_THREADS=8
RUN pip install --no-cache-dir --upgrade pip wheel && \
    pip install --no-cache-dir \
        pandas \
        scipy \
        uncertainties \
        typing \
        pathlib \
        matplotlib \
        numpy && \
        polars>=0.14.1

RUN cd /opt/venv/bin && \
    wget https://temp.aertslab.org/.barcard/create_fragments_file && \
    wget https://temp.aertslab.org/.barcard/coreutils && \
    chmod a+x create_fragments_file coreutils

# install seq (https://github.com/seq-lang/seq/):
ENV SEQ_VERSION=0.11.0
RUN mkdir -p /opt/seq && \
    wget https://github.com/seq-lang/seq/releases/download/v${SEQ_VERSION}/seq-linux-x86_64.tar.gz && \
    tar xzf seq-linux-x86_64.tar.gz --strip-components 1 -C /opt/seq && \
    rm seq-linux-x86_64.tar.gz
ENV PATH="/opt/seq/bin:${PATH}"
ENV OMP_NUM_THREADS=4
ENV SEQ_PYTHON=/usr/lib/x86_64-linux-gnu/libpython3.9.so.1

# install single_cell_toolkit
# https://github.com/aertslab/single_cell_toolkit
RUN git clone --depth=1 https://github.com/aertslab/single_cell_toolkit.git /opt/single_cell_toolkit
ENV seq_root_dir=/opt/seq
ENV PATH="/opt/single_cell_toolkit:/opt/single_cell_toolkit/barcard:${PATH}"

RUN rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/* && \
    ldconfig

