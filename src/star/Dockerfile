FROM debian:buster-slim


ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends apt-utils debconf locales && dpkg-reconfigure locales && \
    # Need to run ps:
    apt-get -y install procps wget less


ARG version="2.7.6a"
RUN cd /usr/local/bin && \
    wget https://github.com/alexdobin/STAR/raw/${version}/bin/Linux_x86_64_static/STAR && \
    chmod +x STAR


