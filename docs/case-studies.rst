Case Studies
=============

Kurmangaliyev Y Z et al., 2019 - Modular transcriptional programs separately define axon and dendrite connectivity
-------------------------------------------------------------------------------------------------------------------

Some links related to the case study:

- Paper: https://elifesciences.org/articles/50822
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126139

Analysis of 10xGenomics Samples
*******************************

BBKNN and SCENIC
++++++++++++++++

The following command was used to generate the config:

.. code:: bash

    nextflow config \
        ~/vib-singlecell-nf/vsn-pipelines \
        -profile sra,cellranger,pcacv,bbknn,dm6,scenic,scenic_use_cistarget_motifs,scenic_use_cistarget_tracks,singularity,qsub \
        > nextflow.config

The generated config is available at the ``vsn-pipelines`` GitHub repository: `examples/kurmangaliyevyz_2019/10x_bbknn_scenic.config`_. You should update the lines commented with " TO EDIT" with the correct information.

.. _`examples/kurmangaliyevyz_2019/10x_bbknn_scenic.config`: https://github.com/vib-singlecell-nf/vsn-pipelines/blob/master/examples/kurmangaliyevyz_2019/10x_bbknn_scenic.config

To start the pipeline, run the following command:

.. code:: bash

    nextflow \
        -C nextflow.config \
        run ~/vib-singlecell-nf/vsn-pipelines \
            -entry sra_cellranger_bbknn_scenic -resume

The resulting loom file is available at `kurmangaliyevyz_2019_10x_bbknn_scenic`_ and is ready to be explored in `SCope <http://scope.aertslab.org/>`_.

.. _`kurmangaliyevyz_2019_10x_bbknn_scenic`: https://cloud.aertslab.org/index.php/s/dpmQyKAW5cWn9RF

Harmony and SCENIC (append mode)
++++++++++++++++++++++++++++++++

.. code:: bash

    nextflow config \
        ~/vib-singlecell-nf/vsn-pipelines \
        -profile tenx,pcacv,harmony,scenic_append_only,singularity \
        > nextflow.config

The generated config is available at the ``vsn-pipelines`` GitHub repository: `examples/kurmangaliyevyz_2019/10x_harmony_scenic.config`_. You should update The lines commented with " TO EDIT" with the correct information.

.. _`examples/kurmangaliyevyz_2019/10x_harmony_scenic.config`: https://github.com/vib-singlecell-nf/vsn-pipelines/blob/master/examples/kurmangaliyevyz_2019/10x_harmony_scenic.config

To start the pipeline, run the following command:

.. code:: bash

    nextflow \
        -C nextflow.config \
        run ~/vib-singlecell-nf/vsn-pipelines \
            -entry harmony_scenic -resume

The resulting loom file is available at `kurmangaliyevyz_2019_harmony_scenic`_ and is ready to be explored in `SCope <http://scope.aertslab.org/>`_.

.. _`kurmangaliyevyz_2019_harmony_scenic`: https://cloud.aertslab.org/index.php/s/92bR4LfLDbtDM8F

Hung R et al., 2019 - A cell atlas of the adult Drosophila midgut
-----------------------------------------------------------------

Some links related to the case study:

- Paper: https://www.pnas.org/content/117/3/1514.abstract
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120537

Analysis of 10xGenomics Samples
*******************************

The following command was used to generate the config:

.. code:: bash

    nextflow config \
        ~/vib-singlecell-nf/vsn-pipelines \
        -profile singularity,sra,cellranger,pcacv,bbknn,scenic \
        > nextflow.config


The generated config is available at the ``vsn-pipelines`` GitHub repository: `examples/hungr_2019/10x_bbknn_scenic.config`_.  You should provide The lines commented with " TO EDIT" with the correct information.

.. _`examples/hungr_2019/10x_bbknn_scenic.config`: https://github.com/vib-singlecell-nf/vsn-pipelines/blob/master/examples/hungr_2019/10x_bbknn_scenic.config

To start the pipeline, run the following command:

.. code:: bash

    nextflow \
        -C nextflow.config \
        run ~/vib-singlecell-nf/vsn-pipelines \
            -entry sra_cellranger_bbknn_scenic


The resulting loom file is available at `hungr_2019_bbknn_scenic.loom`_, and is ready to be explored in `SCope <http://scope.aertslab.org/>`_.

.. _`hungr_2019_bbknn_scenic.loom`: https://cloud.aertslab.org/index.php/s/PeBcfa8ggzbjZRr