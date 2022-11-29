Trim Galore module
==================

This repository contains an implementation of Trim Galore for VIB-SingleCell-NF (VSN) pipelines.
See `here <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`_ for the original source.

To build the Docker image
-------------------------

Image tag format: ```trimgalore-<trimgalore release version>-cutadapt-<cutadapt release version>``.

.. code:: bash

    docker build -t vibsinglecellnf/trimgalore:trimgalore-0.6.7-cutadapt-4.1 .
    podman build -t vibsinglecellnf/trimgalore:trimgalore-0.6.7-cutadapt-4.1 .

This image uses the ``vibsinglecellnf/samtools`` image as a base.

