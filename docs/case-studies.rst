Case Studies
=============

Hung R et al., 2019 - A cell atlas of the adult Drosophila midgut
-----------------------------------------------------------------

Some links related to the case study:

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120537
- Paper: https://www.pnas.org/content/117/3/1514.abstract

Analysis of 10x Samples
************************

The following command was used to generate the config::

    nextflow config \
        ~/vib-singlecell-nf \
        -profile singularity,sra,cellranger,pcacv,bbknn,scenic \
        > nextflow.config


The generated config is available in ``examples/hungr_2019/10x_bbknn_scenic.config``.

To start the pipeline, run the following command::

    nextflow \
    -C nextflow.config \
    run ~/vib-singlecell-nf \
        -entry sra_cellranger_bbknn_scenic


The resulting loom file is available here: ``examples/hungr_2019/10x_bbknn_scenic.loom`` and is ready to be explored in `SCope <http://scope.aertslab.org/>`_.
