FROM qrouchon/r-base-plus:4.0.2
RUN apt-get update

RUN apt-get install -y libgit2-dev libpng-dev && rm -rf /var/lib/apt/lists/*

# Python dependencies
# python3 python3-pip
#RUN pip3 install --upgrade setuptools
#RUN pip3 install numpy scipy scikit-learn numba umap-learn

RUN R -e 'BiocManager::install(c("SingleCellExperiment","multtest"),update = TRUE, ask = FALSE);install.packages(c("Seurat"));devtools::install_github("mojaveazure/seurat-disk");'
