FROM vibsinglecellnf/samtools:1.10
# source: https://github.com/vib-singlecell-nf/bwamaptools/blob/master/Dockerfile.samtools

RUN apk add --virtual build-dependencies git cmake g++

# Install popscle
RUN git clone --depth 1 https://github.com/statgen/popscle.git && \
    cd popscle && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    cp /popscle/bin/popscle /usr/local/bin

# install popscle_helper_tools into this image
# (https://github.com/aertslab/popscle_helper_tools)
RUN git clone --depth 1 https://github.com/aertslab/popscle_helper_tools.git /tmp/popscle_helper_tools && \
    mv /tmp/popscle_helper_tools/*sh /usr/bin

RUN apk del build-dependencies && \
    rm -rf /var/cache/apk/*

