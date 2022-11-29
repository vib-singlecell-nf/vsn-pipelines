FROM vibsinglecellnf/samtools:0.3-1.16.1

# Install popscle.
RUN git clone --depth 1 https://github.com/statgen/popscle.git /tmp/popscle && \
    mkdir -p /tmp/popscle/build && \
    cd /tmp/popscle/build && \
    cmake .. && \
    make && \
    mv /tmp/popscle/bin/popscle /usr/local/bin && \
    rm -rf /tmp/popscle

# Install bedtools.
ENV BEDTOOLS_VERSION 2.30.0
RUN curl -L -o /usr/local/bin/bedtools \
      https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static.binary && \
    chmod a+x /usr/local/bin/bedtools

# Install filter_bam_file_for_popscle_dsc_pileup.sh from popscle_helper_tools into this image.
# (https://github.com/aertslab/popscle_helper_tools)
RUN git clone --depth 1 https://github.com/aertslab/popscle_helper_tools.git /tmp/popscle_helper_tools && \
    mv /tmp/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh /usr/local/bin && \
    rm -rf /tmp/popscle_helper_tools

