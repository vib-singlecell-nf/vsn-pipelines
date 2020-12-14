FROM qrouchon/r-base-plus:4.0.2
LABEL maintainer="QuR <quentin.rouchon@irc.vib-ugent.be>"

RUN apt-get update -qq \
  && apt-get install -y --no-install-recommends libcurl4-openssl-dev libhdf5-dev

RUN R -e 'BiocManager::install(c("scran","DropletUtils","scMerge"),update = TRUE, ask = FALSE)'
