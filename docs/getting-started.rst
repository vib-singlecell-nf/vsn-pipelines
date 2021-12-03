Getting Started
================

Prerequisite
************

Make sure that ``LANG`` and ``LC_ALL`` environment variables have been set. You can use the following command to check this:

.. code:: shell

    locale

If some are not set, you can set them to the default language for instance:

.. code:: shell

    export LANG="C" 
    export LC_ALL="C"

Dependencies
------------
Make sure you have the following software installed,

- Nextflow_
    
    - Currently VSN-Pipelines requires Nextflow version ``21.04.03`` or higher.

- A container system, either of:

    - Docker_
    - Singularity_

**NOTE**: Due to licensing restrictions, to use the cellranger components of VSN you must build and/or provide a container with ``cellranger`` and ``bcl2fastq2`` installed yourself.
A sample ``Dockerfile`` can be found in ``./src/cellranger/``, you must download ``bcl2fastq2`` from the Illumina_ website, and cellranger from the `10x Genomics`_ website yourself to build this container.

.. _Nextflow: https://www.nextflow.io/
.. _Docker: https://docs.docker.com/
.. _Singularity: https://www.sylabs.io/singularity/
.. _Illumina: https://emea.support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html
.. _`10x Genomics`: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger

Quick start
***********

To run a quick test of the single sample analysis pipeline, we can use the 1k PBMC datasets provided by 10x Genomics.
This will take only **~3min** to run.

1. The data first needs to be downloaded (instructions can be found `here <https://github.com/vib-singlecell-nf/vsn-pipelines/tree/master/data>`_).

2. Next, update to the latest pipeline version::

    nextflow pull vib-singlecell-nf/vsn-pipelines

3. Next, generate a config file using the standard settings for the test data, and the appropriate profiles (e.g., replace ``singularity`` with ``docker`` if necessary)::

    nextflow config vib-singlecell-nf/vsn-pipelines \
        -profile tenx,singularity,single_sample > single_sample.config

4. The test pipeline can now be run using the config file just generated, specifying the ``single_sample`` workflow as an entrypoint::

    nextflow -C single_sample.config \
        run vib-singlecell-nf/vsn-pipelines \
            -entry single_sample

Example Output
--------------

.. code:: shell

    $ nextflow -C nextflow_new_structure.config run $VSN -entry single_sample
    N E X T F L O W  ~  version 21.04.3
    Launching `/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/GitHub/vib-singlecell-nf/vsn-pipelines/main.nf` [loving_shockley] - revision: ba1dedbf51
    executor >  local (75)
    [1d/3b5a55] process > single_sample:SINGLE_SAMPLE:SC__FILE_CONVERTER (2)                                                                                [100%] 2 of 2 _
    [2d/2152cf] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:QC_FILTER:SC__SCANPY__COMPUTE_QC_STATS (2)                                      [100%] 2 of 2 _
    [48/bce024] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:QC_FILTER:SC__SCANPY__CELL_FILTER (2)                                           [100%] 2 of 2 _
    [60/d42cdf] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:QC_FILTER:SC__SCANPY__GENE_FILTER (2)                                           [100%] 2 of 2 _
    [4b/bb2635] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:QC_FILTER:GENERATE_DUAL_INPUT_REPORT:SC__SCANPY__GENERATE_DUAL_INPUT_REPORT (2) [100%] 2 of 2 _
    [64/add548] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:QC_FILTER:GENERATE_DUAL_INPUT_REPORT:SC__SCANPY__REPORT_TO_HTML (2)             [100%] 2 of 2 _
    [c6/4d8a66] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:NORMALIZE_TRANSFORM:SC__SCANPY__NORMALIZATION (2)                               [100%] 2 of 2 _
    [01/8ba0d2] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:NORMALIZE_TRANSFORM:PUBLISH_H5AD_NORMALIZED:COMPRESS_HDF5 (2)                   [100%] 2 of 2 _
    [b3/ec4712] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:NORMALIZE_TRANSFORM:PUBLISH_H5AD_NORMALIZED:SC__PUBLISH (2)                     [100%] 2 of 2 _
    [1e/35bb2e] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:NORMALIZE_TRANSFORM:SC__SCANPY__DATA_TRANSFORMATION (2)                         [100%] 2 of 2 _
    [14/adfd7a] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:HVG_SELECTION:SC__SCANPY__FIND_HIGHLY_VARIABLE_GENES (2)                        [100%] 2 of 2 _
    [01/9c8a26] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:HVG_SELECTION:SC__SCANPY__SUBSET_HIGHLY_VARIABLE_GENES (2)                      [100%] 2 of 2 _
    [dc/027334] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:HVG_SELECTION:SC__SCANPY__FEATURE_SCALING (2)                                   [100%] 2 of 2 _
    [8d/05ce2f] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:HVG_SELECTION:PUBLISH_H5AD_HVG_SCALED:COMPRESS_HDF5 (2)                         [100%] 2 of 2 _
    [0b/6d50b0] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:HVG_SELECTION:PUBLISH_H5AD_HVG_SCALED:SC__PUBLISH (2)                           [100%] 2 of 2 _
    [c1/f799be] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:HVG_SELECTION:GENERATE_REPORT:SC__SCANPY__GENERATE_REPORT (2)                   [100%] 2 of 2 _
    [c9/ae0cd9] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:HVG_SELECTION:GENERATE_REPORT:SC__SCANPY__REPORT_TO_HTML (2)                    [100%] 2 of 2 _
    [a2/0a7824] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:DIM_REDUCTION_PCA:SC__SCANPY__DIM_REDUCTION__PCA (2)                            [100%] 2 of 2 _
    [d6/1407b3] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:NEIGHBORHOOD_GRAPH:SC__SCANPY__NEIGHBORHOOD_GRAPH (2)                           [100%] 2 of 2 _
    [b7/8ab962] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:DIM_REDUCTION_TSNE_UMAP:SC__SCANPY__DIM_REDUCTION__TSNE (2)                     [100%] 2 of 2 _
    [ee/485413] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:DIM_REDUCTION_TSNE_UMAP:SC__SCANPY__DIM_REDUCTION__UMAP (2)                     [100%] 2 of 2 _
    [ba/2bfd23] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:DIM_REDUCTION_TSNE_UMAP:GENERATE_REPORT:SC__SCANPY__GENERATE_REPORT (2)         [100%] 2 of 2 _
    [10/a429ce] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:DIM_REDUCTION_TSNE_UMAP:GENERATE_REPORT:SC__SCANPY__REPORT_TO_HTML (2)          [100%] 2 of 2 _
    [06/3412cd] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:SC__SCANPY__CLUSTERING (2)                               [100%] 2 of 2 _
    [23/3d82c4] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:GENERATE_REPORT:SC__SCANPY__GENERATE_REPORT (2)          [100%] 2 of 2 _
    [bb/c9e11f] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:GENERATE_REPORT:SC__SCANPY__REPORT_TO_HTML (2)           [100%] 2 of 2 _
    [1c/2026be] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:SC__SCANPY__MARKER_GENES (2)                             [100%] 2 of 2 _
    [57/13f0a8] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:UTILS__GENERATE_WORKFLOW_CONFIG_REPORT                                          [100%] 1 of 1 _
    [60/8a3231] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:SC__SCANPY__MERGE_REPORTS (2)                                                   [100%] 2 of 2 _
    [cb/de1a4d] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:SC__SCANPY__REPORT_TO_HTML (2)                                                  [100%] 2 of 2 _
    [3f/265503] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:FINALIZE:SC__H5AD_TO_FILTERED_LOOM (2)                                          [100%] 2 of 2 _
    [1f/de67e8] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:FINALIZE:FILE_CONVERTER_TO_SCOPE:SC__H5AD_TO_LOOM (2)                           [100%] 2 of 2 _
    [2a/10d5a2] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:FINALIZE:FILE_CONVERTER_TO_SCANPY:SC__H5AD_MERGE (2)                            [100%] 2 of 2 _
    [35/ce7256] process > single_sample:SINGLE_SAMPLE:SCANPY__SINGLE_SAMPLE:PUBLISH:SC__PUBLISH (2)                                                         [100%] 2 of 2 _
    [6f/1df294] process > single_sample:PUBLISH_SINGLE_SAMPLE_SCOPE:COMPRESS_HDF5 (2)                                                                       [100%] 2 of 2 _
    [90/a9e563] process > single_sample:PUBLISH_SINGLE_SAMPLE_SCOPE:SC__PUBLISH (2)                                                                         [100%] 2 of 2 _
    [23/d62b2e] process > single_sample:PUBLISH_SINGLE_SAMPLE_SCANPY:COMPRESS_HDF5 (2)                                                                      [100%] 2 of 2 _
    [d0/5834be] process > single_sample:PUBLISH_SINGLE_SAMPLE_SCANPY:SC__PUBLISH (1)                                                                        [100%] 2 of 2 _
    WARN: To render the execution DAG in the required format it is required to install Graphviz -- See http://www.graphviz.org for more info.
    Completed at: 26-Aug-2021 15:41:37
    Duration    : 2m 29s
    CPU hours   : 0.4
    Succeeded   : 75


Output
------

The pipelines will generate 3 types of results in the output directory (`params.global.outdir`), by default ``out/``

- ``data``: contains the workflow output file (in h5ad format), plus symlinks to all the intermediate files.
- ``loom``: contains final loom files which can be imported inside SCope visualization tool for further visualization of the results.
- ``notebooks``: contains all the notebooks generated along the pipeline (e.g.: Quality control report)
- ``pipeline_reports``: Nextflow dag, execution, timeline, and trace reports

For a full list of the pipelines available please see the `pipelines <./pipelines.html>`_ page.


----


Further pipeline configuration details
**************************************

This pipeline can be fully configured and run on custom data with a few steps.
The recommended method is to first run ``nextflow config ...`` to generate a complete config file (with the default parameters) in your working directory.
The tool-specific parameters, as well as Docker/Singularity profiles, are included when specifying the appropriate profiles to ``nextflow config``.

1. First, update to the latest pipeline version (this will update the Nextflow cache of the repository, typically located in ``~/.nextflow/assets/vib-singlecell-nf/``)::

    nextflow pull vib-singlecell-nf/vsn-pipelines


2. Next, a config file needs to be generated.
   This step will merge parameters from multiple profiles together to create a master config which specifies **all** parameters used by the pipeline.
   In this example, these are ``tenx`` for the input data, ``singularity`` to use the Singularity system (replace with ``docker`` if necessary), and ``single_sample`` to load the defaults for the single sample pipeline.
   In your working directory, run ``nextflow config ...`` with the appropriate profiles::

    nextflow config vib-singlecell-nf/vsn-pipelines \
        -profile tenx,singularity,single_sample > single_sample.config



3. Now, edits can be made to ``single_sample.config``.
   Generally, the default values are acceptable to use for a first pass, but certain variables (input directory, etc.) need to be changed.

   In particular, the following parameters are frequently modified in practice:

    * ``params.global.project_name``: a project name which will be included in some of the output file names.
    * ``params.data.tenx.cellranger_mex``, which should point to the ``outs/`` folder generated by Cell Ranger (if using 10x data). See :ref:`Information on using 10x Genomics datasets<using_10x_datasets>` for additional info.
    * Filtering parameters (``params.tools.scanpy.filter``): filtering parameters, which will be applied to all samples, can be set here: min/max genes, mitochondrial read fraction, and min cells. See :ref:`Multi-sample parameters` for additional info on how to specify sample-specific parameters.
    * Louvain cluster resolution: ``params.tools.scanpy.clustering.resolution``.
    * :ref:`Cell-<Cell annotation>` and :ref:`sample-<Sample annotation>` level annotations are also possible.


4. Run the workflow using the new config file (using ``-C`` is recommended to use **only** this file), specifying the proper workflow as the entry point::

    nextflow -C single_sample.config \
        run vib-singlecell-nf/vsn-pipelines \
        -entry single_sample


Additional resources for running on custom data
-----------------------------------------------

- `Input file formats available <./input_formats.html>`_.
- `Available pipelines <./pipelines.html>`_.
- `Advanced features for customizing pipelines <./features.html>`_.

Finally, see the list of case studies with specific examples and full config files at `VSN-Pipelines-examples <https://vsn-pipelines-examples.readthedocs.io/en/latest/>`_.

