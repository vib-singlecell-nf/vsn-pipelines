FROM vibsinglecellnf/samtools:0.2-1.12

RUN echo 'deb http://deb.debian.org/debian testing main' > /etc/apt/sources.list.d/testing.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        cmake

# Install popscle
RUN git clone --depth 1 https://github.com/statgen/popscle.git /tmp/popscle && \
    mkdir -p /tmp/popscle/build && \
    cd /tmp/popscle/build && \
    cmake .. && \
    make && \
    cp /tmp/popscle/bin/popscle /usr/local/bin

# install bedtools
ENV BEDTOOLS_VERSION 2.30.0
RUN curl -L -o /usr/local/bin/bedtools \
    https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static.binary && \
    chmod a+x /usr/local/bin/bedtools

# install popscle_helper_tools into this image
# (https://github.com/aertslab/popscle_helper_tools)
RUN git clone --depth 1 https://github.com/aertslab/popscle_helper_tools.git /tmp/popscle_helper_tools && \
    mv /tmp/popscle_helper_tools/*sh /usr/local/bin

