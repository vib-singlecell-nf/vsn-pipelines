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

A repository of pipelines for single-cell data in Nextflow DSL2.

Full documentation available on `Read the Docs <https://vsn-pipelines.readthedocs.io/en/latest/>`_

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
