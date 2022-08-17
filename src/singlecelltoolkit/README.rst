
single_cell_toolkit template
============================

This repository contains an implementation of single_cell_toolkit for VIB-SingleCell-NF (VSN) pipelines.
See `aertslab/single_cell_toolkit <https://github.com/aertslab/single_cell_toolkit>`_ for the original source.

To build the Docker image
-------------------------

Image tag format: ``<date of latest git commit>-<short hash of latest git commit>``.

.. code:: bash

    docker build -t vibsinglecellnf/singlecelltoolkit:2022-07-07-0638c1d .
    podman build -t vibsinglecellnf/singlecelltoolkit:2022-07-07-0638c1d .

This image uses the ``vibsinglecellnf/samtools`` image as a base.

