FROM vibsinglecellnf/samtools:1.11
# source: https://github.com/vib-singlecell-nf/bwamaptools/blob/master/Dockerfile.samtools

RUN apk add --virtual build-dependencies git

# install single_cell_toolkit
# https://github.com/aertslab/single_cell_toolkit
RUN git clone https://github.com/aertslab/single_cell_toolkit.git /tmp/single_cell_toolkit && \
    mv /tmp/single_cell_toolkit/*sh /usr/bin

RUN apk del build-dependencies && \
    rm -rf /var/cache/apk/*

