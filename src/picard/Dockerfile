FROM adoptopenjdk/openjdk8:jre8u-nightly

ARG build_command=shadowJar
ENV jar_name=picard.jar
ENV picard_version=2.21.1

# Install ant, git for building
RUN apt-get update && \
    apt-get install wget && \
    apt-get clean autoclean && \
    apt-get autoremove -y

WORKDIR /usr/local/bin
RUN wget https://github.com/broadinstitute/picard/releases/download/${picard_version}/${jar_name}
