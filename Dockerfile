FROM rocker/r-ver:3.6.0

RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev libxml2-dev zlib1g-dev libhdf5-dev && \
    apt-get install -y libssl-dev && \
    R -e "install.packages('devtools')" && \
    R -e "install.packages('RCurl')" && \
    R -e "install.packages('rvest')" && \
    R -e "devtools::install_github('dweemx/flybaseR')"

RUN cd / && \
   rm -rf /tmp/* && \
   apt-get autoremove -y && \
   apt-get autoclean -y && \
   rm -rf /var/cache/apt/* && \
   rm -rf /var/lib/apt/lists/* && \
   apt-get clean