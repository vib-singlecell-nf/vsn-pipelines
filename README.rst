VSN-Pipelines
==============

|Nextflow| |Gitter| |ReadTheDocs|

.. |ReadTheDocs| image:: https://readthedocs.org/projects/vsn-pipelines/badge/?version=latest
    :target: https://vsn-pipelines.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |Nextflow| image:: https://img.shields.io/badge/nextflow-19.12.0-brightgreen.svg
    :target: https://www.nextflow.io/
    :alt: Nextflow

.. |Gitter| image:: https://badges.gitter.im/vib-singlecell-nf/community.svg
    :target: https://gitter.im/vib-singlecell-nf/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
    :alt: Gitter

.. |Zenodo| image:: https://zenodo.org/badge/199477571.svg
    :target: https://zenodo.org/badge/latestdoi/199477571
    :alt: Zenodo https://zenodo.org/badge/latestdoi/199477571


|single_sample| |single_sample_scenic| |scenic| |scenic_multiruns| |single_sample_scenic_multiruns| |bbknn| |bbknn_scenic| |harmony| |mnncorrect|

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

.. |bbknn| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/bbknn/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#bbknn-bbknn
    :alt: BBKNN Pipeline

.. |bbknn_scenic| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/bbknn_scenic/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#bbknn-scenic
    :alt: BBKNN SCENIC Pipeline

.. |harmony| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/harmony/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#harmony-harmony
    :alt: Harmony Pipeline

.. |mnncorrect| image:: https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/mnncorrect/badge.svg
    :target: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#mnncorrect-mnncorrect
    :alt: MNN-correct Pipeline

A repository of pipelines for single-cell data in Nextflow DSL2.

A quick tour of the VSN pipelines ? Please read `Quick Start <https://vsn-pipelines.readthedocs.io/en/latest/getting-started.html#quick-start>`_.

Full documentation available on `Read the Docs <https://vsn-pipelines.readthedocs.io/en/latest/>`_

If VSN-Pipelines is useful for your research, consider citing:

- VSN-Pipelines v0.14.0: `10.5281/zenodo.3703109 <https://doi.org/10.5281/zenodo.3703109>`_.
- VSN-Pipelines All Versions (will always resolve to the latest one): `10.5281/zenodo.3703108 <https://doi.org/10.5281/zenodo.3703108>`_.

This main repo contains multiple workflows for analyzing single cell transcriptomics data, and depends on a number of tools, which are organized into submodules within the VIB-Singlecell-NF_ organization.
Currently available workflows include:

.. _VIB-Singlecell-NF: https://github.com/vib-singlecell-nf

- **Cell Ranger**: processes 10x Chromium data to align reads to generate an expression counts matrix.
- **DropSeq**: processes Drop-seq data from read alignment to expression counts.
- **Single sample workflows**: perform a "best practices" scRNA-seq analysis. Multiple samples can be run in parallel, treating each sample separately.
- **Multi-sample workflows**: perform a "best practices" scRNA-seq analysis on a merged and batch-corrected group of samples. Available batch correction methods include:

    - **BBKNN**
    - **mnnCorrect**
    - **Harmony**

* **GRN inference**:

    * The pySCENIC_ implementation of the SCENIC_ workflow is integrated here and can be run in conjunction with any of the above workflows.

.. _pySCENIC: https://github.com/aertslab/pySCENIC
.. _SCENIC: https://aertslab.org/#scenic

The output of each of the main workflows is a loom_-format file, which is ready for import into the interactive single-cell web visualization tool SCope_.
In addition, data is also output in h5ad format, and reports are generated for the major pipeline steps.

.. _loom: http://loompy.org/
.. _SCope: http://scope.aertslab.org/
