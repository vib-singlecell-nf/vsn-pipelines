FROM ncbi/sra-toolkit AS compile-image

# set the environment variables
ENV pigz_version 2.4

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
RUN wget https://github.com/madler/pigz/archive/v${pigz_version}.zip && \
   unzip v${pigz_version}.zip && \
   cd pigz-${pigz_version} && \
   make

FROM ncbi/sra-toolkit AS build-image
ENV pigz_version 2.4

COPY --from=compile-image /usr/local/bin/pigz-${pigz_version}/pigz /usr/local/bin
COPY --from=compile-image /usr/local/bin/pigz-${pigz_version}/unpigz /usr/local/bin
