VSN-Pipelines
==============

A repository of pipelines for single-cell data analysis in Nextflow DSL2.

|VSN-Pipelines| |ReadTheDocs| |Zenodo| |Gitter| |Nextflow|


**Full documentation** is available on `Read the Docs <https://vsn-pipelines.readthedocs.io/en/latest/>`_, or take a look at the `Quick Start <https://vsn-pipelines.readthedocs.io/en/latest/getting-started.html#quick-start>`_ guide.

This main repo contains multiple workflows for analyzing single cell transcriptomics data, and depends on a number of tools, which are organized into subfolders within the ``src/`` directory.
The VIB-Singlecell-NF_ organization contains this main repo along with a collection of example runs (`VSN-Pipelines-examples <https://vsn-pipelines-examples.readthedocs.io/en/latest/>`_).
Currently available workflows are listed below.

If VSN-Pipelines is useful for your research, consider citing:

- VSN-Pipelines All Versions (latest): `10.5281/zenodo.3703108 <https://doi.org/10.5281/zenodo.3703108>`_.

Raw Data Processing Workflows
-----------------------------

These are set up to run Cell Ranger and DropSeq pipelines.

.. list-table:: Raw Data Processing Workflows
    :widths: 15 10 30
    :header-rows: 1

    * - Pipeline / Entrypoint
      - Purpose
      - Documentation
    * - cellranger
      - Process 10x Chromium data
      - cellranger_
    * - demuxlet_freemuxlet
      - Demultiplexing
      - demuxlet_freemuxlet_
    * - nemesh
      - Process Drop-seq data
      - nemesh_

.. _cellranger: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#cellranger
.. _demuxlet_freemuxlet: https://vsn-pipelines.readthedocs.io/en/develop/pipelines.html#demuxlet-freemuxlet
.. _nemesh: https://vsn-pipelines.readthedocs.io/en/develop/pipelines.html#nemesh


Single Sample Workflows
-----------------------

The **Single Sample Workflows** perform a "best practices" scRNA-seq analysis. Multiple samples can be run in parallel, treating each sample separately.

.. list-table:: Single Sample Workflows
    :header-rows: 1

    * - Pipeline / Entrypoint
      - Purpose
      - Documentation
    * - single_sample
      - Independent samples
      - |single_sample|
    * - single_sample_scenic
      - Ind. samples + SCENIC
      - |single_sample_scenic|
    * - scenic
      - SCENIC GRN inference
      - |scenic|
    * - scenic_multiruns
      - SCENIC run multiple times
      - |scenic_multiruns|
    * - single_sample_scenic_multiruns
      - Ind. samples + multi-SCENIC
      - |single_sample_scenic_multiruns|
    * - single_sample_scrublet
      - Ind. samples + Scrublet
      - |single_sample_scrublet|
    * - decontx
      - DecontX
      - |decontx|
    * - single_sample_decontx
      - Ind. samples + DecontX
      - |single_sample_decontx|
    * - single_sample_decontx_scrublet
      - Ind. samples + DecontX + Scrublet
      - |single_sample_decontx_scrublet|


Sample Aggregation Workflows
----------------------------

**Sample Aggregation Workflows**: perform a "best practices" scRNA-seq analysis on a merged and batch-corrected group of samples. Available batch correction methods include BBKNN, mnnCorrect, and Harmony.

.. list-table:: Sample Aggregation Pipelines
    :widths: 15 10 30
    :header-rows: 1

    * - Pipeline / Entrypoint
      - Purpose
      - Documentation
    * - bbknn
      - Sample aggregation + BBKNN
      - |bbknn|
    * - bbknn_scenic
      - BBKNN + SCENIC
      - |bbknn_scenic|
    * - harmony
      - Sample aggregation + Harmony
      - |harmony|
    * - harmony_scenic
      - Harmony + SCENIC
      - |harmony_scenic|
    * - mnncorrect
      - Sample aggregation + mnnCorrect
      - |mnncorrect|


----

In addition, the pySCENIC_ implementation of the SCENIC_ workflow is integrated here and can be run in conjunction with any of the above workflows.
The output of each of the main workflows is a loom_-format file, which is ready for import into the interactive single-cell web visualization tool SCope_.
In addition, data is also output in h5ad format, and reports are generated for the major pipeline steps.

scATAC-seq workflows
--------------------

Single cell ATAC-seq processing steps are now included in VSN Pipelines.
Currently, a preprocesing workflow is available, which will take fastq inputs, apply barcode correction, read trimming, bwa mapping, and output bam and fragments files for further downstream analysis.
See `here <https://vsn-pipelines.readthedocs.io/en/latest/scatac-seq.html>`_ for complete documentation.


.. |VSN-Pipelines| image:: https://img.shields.io/github/v/release/vib-singlecell-nf/vsn-pipelines
    :target: https://github.com/vib-singlecell-nf/vsn-pipelines/releases
    :alt: GitHub release (latest by date)

.. |ReadTheDocs| image:: https://readthedocs.org/projects/vsn-pipelines/badge/?version=latest
    :target: https://vsn-pipelines.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |Nextflow| image:: https://img.shields.io/badge/nextflow-21.04.3-brightgreen.svg
    :target: https://www.nextflow.io/
    :alt: Nextflow

.. |Gitter| image:: https://badges.gitter.im/vib-singlecell-nf/community.svg
    :target: https://gitter.im/vib-singlecell-nf/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
    :alt: Gitter

.. |Zenodo| image:: https://zenodo.org/badge/199477571.svg
    :target: https://zenodo.org/badge/latestdoi/199477571
    :alt: Zenodo

.. _VIB-Singlecell-NF: https://github.com/vib-singlecell-nf
.. _pySCENIC: https://github.com/aertslab/pySCENIC
.. _SCENIC: https://aertslab.org/#scenic
.. _loom: http://loompy.org/
.. _SCope: http://scope.aertslab.org/

.. |single_sample| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/single_sample/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#single-sample-single-sample
    :alt: Single-sample Pipeline

.. |single_sample_scenic| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/single_sample_scenic/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#single-sample-scenic-single-sample-scenic
    :alt: Single-sample SCENIC Pipeline

.. |scenic| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/scenic/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#scenic-scenic
    :alt: SCENIC Pipeline

.. |scenic_multiruns| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/scenic_multiruns/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#scenic-multiruns-scenic-multiruns-single-sample-scenic-multiruns
    :alt: SCENIC Multi-runs Pipeline

.. |single_sample_scenic_multiruns| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/single_sample_scenic_multiruns/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#scenic-multiruns-scenic-multiruns-single-sample-scenic-multiruns
    :alt: Single-sample SCENIC Multi-runs Pipeline

.. |single_sample_scrublet| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/single_sample_scrublet/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#single-sample-scrublet-single-sample-scrublet
    :alt: Single-sample Scrublet Pipeline

.. |decontx| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/decontx/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#decontx-decontx
    :alt: DecontX Pipeline

.. |single_sample_decontx| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/single_sample_decontx/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#single-sample-decontx-single-sample-decontx
    :alt: Single-sample DecontX Pipeline

.. |single_sample_decontx_scrublet| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/single_sample_decontx_scrublet/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#single-sample-decontx-scrublet-single-sample-decontx-scrublet
    :alt: Single-sample DecontX Scrublet Pipeline

.. |bbknn| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/bbknn/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#bbknn-bbknn
    :alt: BBKNN Pipeline

.. |bbknn_scenic| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/bbknn_scenic/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#bbknn-scenic
    :alt: BBKNN SCENIC Pipeline

.. |harmony| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/harmony/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#harmony-harmony
    :alt: Harmony Pipeline

.. |harmony_scenic| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/harmony_scenic/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#harmony-scenic
    :alt: Harmony SCENIC Pipeline

.. |mnncorrect| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/mnncorrect/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#mnncorrect-mnncorrect
    :alt: MNN-correct Pipeline

