
VSN-Pipelines popscle
======================

This is a repository for the popslce module of the VIB-SingleCell-NF (VSN) pipelines.

Current Status
---------------

This module currently has two workflows ``freemuxlet`` and ``demuxlet``. 
Both of these workflows expect an input channel consisting of a tuple where
element 1 is the sampleID and element 2 is the output folder of a 10X run.

Currently the workflows are fixed to the filtered matrices.

To build the Docker image
-------------------------

Image tag format: ``<date of latest git commit>-<short hash of latest git commit>``.

.. code:: bash

    docker build -t vibsinglecellnf/popscle:2021-05-05-da70fc7 .

This image uses the ``vibsinglecellnf/samtools`` image as a base.

Acknowledgements
----------------

This module implements functionality developed by Gert Hulselmens designed to
speed up the running time of dsc-pileup. The `filter_bam_file_for_popscle_dsc_pileup`_
script can lead to speedups of 5-10x depending on the input data.

.. _`filter_bam_file_for_popscle_dsc_pileup`: https://github.com/aertslab/popscle_helper_tools

