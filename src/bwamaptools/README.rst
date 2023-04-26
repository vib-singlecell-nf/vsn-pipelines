
BWA maptools module
===================

This repository contains an implementation of BWA for VIB-SingleCell-NF (VSN) pipelines, along with several supporing tools (htslib, samtools).
See `lh3/bwa <https://github.com/lh3/bwa>`_ for the original source.

To build the Docker image
-------------------------

Image tag format: ``<bwa-mem2>-<bwa-mem2 release version>``.

.. code:: bash

    docker build -t vibsinglecellnf/bwamaptools:bwa-mem2-2.2.1-zlibng-2.0.6 .
    podman build -t vibsinglecellnf/bwamaptools:bwa-mem2-2.2.1-zlibng-2.0.6 .

This image uses the ``vibsinglecellnf/samtools`` image as a base.

