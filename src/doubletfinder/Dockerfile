FROM qrouchon/r-base-plus:4.0.2
RUN apt-get update -qq

RUN R -e 'install.packages(c("fields","Seurat,tidyverse"))'
RUN R -e 'devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")'
