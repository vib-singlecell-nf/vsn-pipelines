FROM ubuntu:18.04 AS compile-image

# set the environment variables
ENV fastp_version 0.20.0

RUN apt-get update -y && apt-get install -y \
   build-essential \
   libnss-sss \
   curl \
   wget \
   unzip \
   zlib1g-dev && \
   rm -rf /var/cache/apt/* && \
   rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local/bin
RUN wget https://github.com/OpenGene/fastp/archive/v${fastp_version}.zip && \
   unzip v${fastp_version}.zip && \
   cd fastp-${fastp_version} && \
   make && make install

FROM ubuntu:18.04 AS build-image
RUN apt-get -y update && \
    apt-get -y install zlib1g && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

COPY --from=compile-image /usr/local/bin/fastp /usr/local/bin
