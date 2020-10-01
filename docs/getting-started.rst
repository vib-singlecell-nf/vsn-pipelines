Getting Started
================

Prerequisite
------------

Make sure that ``LANG`` and ``LC_ALL`` environment variables have been set. You can use the following command to check this:

.. code:: shell

    locale

If some are not set, you can set them to the default language for instance:

.. code:: shell

    export LANG="C" 
    export LC_ALL="C"

Dependencies
^^^^^^^^^^^^
Make sure you have the following software installed,

- Nextflow_
- A container system, either of:

    - Docker_
    - Singularity_

.. _Nextflow: https://www.nextflow.io/
.. _Docker: https://docs.docker.com/
.. _Singularity: https://www.sylabs.io/singularity/

Quick start
-----------

To run a quick test of the single sample analysis pipeline, we can use the 1k PBMC datasets provided by 10x Genomics.
This will take only **~3min** to run.

1. The data first needs to be downloaded (instructions can be found here_)

.. _here: ../data/README.md

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
^^^^^^^^^^^^^^

.. code:: shell

    $ nextflow -C single_sample.config run vib-singlecell-nf/vsn-pipelines -entry single_sample

    N E X T F L O W  ~  version 19.12.0-edge
    Launching `/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/GitHub/vib-singlecell-nf/main.nf` [nice_engelbart] - revision: 0096df9054
    WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE
    executor >  local (59)
    [0c/d33a4e] process > single_sample:SINGLE_SAMPLE:UTILS__GENERATE_WORKFLOW_CONFIG_REPORT                                          [100%] 1 of 1 ✔
    [17/ab2b39] process > single_sample:SINGLE_SAMPLE:QC_FILTER:SC__FILE_CONVERTER (1)                                                [100%] 2 of 2 ✔
    [e4/84f688] process > single_sample:SINGLE_SAMPLE:QC_FILTER:SC__SCANPY__COMPUTE_QC_STATS (2)                                      [100%] 2 of 2 ✔
    [1b/daa1c3] process > single_sample:SINGLE_SAMPLE:QC_FILTER:SC__SCANPY__GENE_FILTER (2)                                           [100%] 2 of 2 ✔
    [fc/8653d0] process > single_sample:SINGLE_SAMPLE:QC_FILTER:SC__SCANPY__CELL_FILTER (2)                                           [100%] 2 of 2 ✔
    [9d/ebeff9] process > single_sample:SINGLE_SAMPLE:QC_FILTER:GENERATE_DUAL_INPUT_REPORT:SC__SCANPY__GENERATE_DUAL_INPUT_REPORT (2) [100%] 2 of 2 ✔
    [87/e13dd0] process > single_sample:SINGLE_SAMPLE:QC_FILTER:GENERATE_DUAL_INPUT_REPORT:SC__SCANPY__REPORT_TO_HTML (2)             [100%] 2 of 2 ✔
    [a6/867a4a] process > single_sample:SINGLE_SAMPLE:NORMALIZE_TRANSFORM:SC__SCANPY__NORMALIZATION (2)                               [100%] 2 of 2 ✔
    [07/8e63b1] process > single_sample:SINGLE_SAMPLE:NORMALIZE_TRANSFORM:SC__SCANPY__DATA_TRANSFORMATION (2)                         [100%] 2 of 2 ✔
    [c1/07c18c] process > single_sample:SINGLE_SAMPLE:HVG_SELECTION:SC__SCANPY__FIND_HIGHLY_VARIABLE_GENES (2)                        [100%] 2 of 2 ✔
    [e9/53e204] process > single_sample:SINGLE_SAMPLE:HVG_SELECTION:SC__SCANPY__SUBSET_HIGHLY_VARIABLE_GENES (2)                      [100%] 2 of 2 ✔
    [0b/e7ae8c] process > single_sample:SINGLE_SAMPLE:HVG_SELECTION:SC__SCANPY__FEATURE_SCALING (2)                                   [100%] 2 of 2 ✔
    [5d/52236c] process > single_sample:SINGLE_SAMPLE:HVG_SELECTION:GENERATE_REPORT:SC__SCANPY__GENERATE_REPORT (2)                   [100%] 2 of 2 ✔
    [71/5d6559] process > single_sample:SINGLE_SAMPLE:HVG_SELECTION:GENERATE_REPORT:SC__SCANPY__REPORT_TO_HTML (2)                    [100%] 2 of 2 ✔
    [8c/1b4cc9] process > single_sample:SINGLE_SAMPLE:DIM_REDUCTION_PCA:SC__SCANPY__DIM_REDUCTION__PCA (2)                            [100%] 2 of 2 ✔
    [7b/d423f7] process > single_sample:SINGLE_SAMPLE:NEIGHBORHOOD_GRAPH:SC__SCANPY__NEIGHBORHOOD_GRAPH (2)                           [100%] 2 of 2 ✔
    [9b/3a10d2] process > single_sample:SINGLE_SAMPLE:DIM_REDUCTION_TSNE_UMAP:SC__SCANPY__DIM_REDUCTION__TSNE (2)                     [100%] 2 of 2 ✔
    [5f/2c6325] process > single_sample:SINGLE_SAMPLE:DIM_REDUCTION_TSNE_UMAP:SC__SCANPY__DIM_REDUCTION__UMAP (2)                     [100%] 2 of 2 ✔
    [ff/b5c6ef] process > single_sample:SINGLE_SAMPLE:DIM_REDUCTION_TSNE_UMAP:GENERATE_REPORT:SC__SCANPY__GENERATE_REPORT (2)         [100%] 2 of 2 ✔
    [b6/86bc36] process > single_sample:SINGLE_SAMPLE:DIM_REDUCTION_TSNE_UMAP:GENERATE_REPORT:SC__SCANPY__REPORT_TO_HTML (2)          [100%] 2 of 2 ✔
    [1a/2fec91] process > single_sample:SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:SC__SCANPY__CLUSTERING (2)                               [100%] 2 of 2 ✔
    [38/8a814b] process > single_sample:SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:GENERATE_REPORT:SC__SCANPY__GENERATE_REPORT (2)          [100%] 2 of 2 ✔
    [35/530dcf] process > single_sample:SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:GENERATE_REPORT:SC__SCANPY__REPORT_TO_HTML (2)           [100%] 2 of 2 ✔
    [05/3e201e] process > single_sample:SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:SC__SCANPY__MARKER_GENES (2)                             [100%] 2 of 2 ✔
    [04/ad44c6] process > single_sample:SINGLE_SAMPLE:SC__H5AD_TO_FILTERED_LOOM (2)                                                   [100%] 2 of 2 ✔
    [46/47cac6] process > single_sample:SINGLE_SAMPLE:FILE_CONVERTER:SC__H5AD_TO_LOOM (2)                                             [100%] 2 of 2 ✔
    [33/640ffa] process > single_sample:SINGLE_SAMPLE:FILE_CONVERTER:COMPRESS_HDF5 (2)                                                [100%] 2 of 2 ✔
    [77/87b596] process > single_sample:SINGLE_SAMPLE:PUBLISH (2)                                                            [100%] 2 of 2 ✔
    [61/82bf98] process > single_sample:SINGLE_SAMPLE:SC__SCANPY__MERGE_REPORTS (1)                                                   [100%] 2 of 2 ✔
    [5a/26ce75] process > single_sample:SINGLE_SAMPLE:SC__SCANPY__REPORT_TO_HTML (2)                                                  [100%] 2 of 2 ✔

    ------------------------------------------------------------------
    Converting 1k_pbmc_v2_chemistry.SC__SCANPY__MARKER_GENES.h5ad to 1k_pbmc_v2_chemistry.SC__SCANPY__MARKER_GENES.loom
    (w/ additional compression)...
    ------------------------------------------------------------------


    ------------------------------------------------------------------
    Converting 1k_pbmc_v3_chemistry.SC__SCANPY__MARKER_GENES.h5ad to 1k_pbmc_v3_chemistry.SC__SCANPY__MARKER_GENES.loom
    (w/ additional compression)...
    ------------------------------------------------------------------

    WARN: To render the execution DAG in the required format it is required to install Graphviz -- See http://www.graphviz.org for more info.
    Completed at: 25-Feb-2020 12:31:44
    Duration    : 2m 15s
    CPU hours   : 0.1
    Succeeded   : 59


The pipelines will generate 3 types of results in the output directory (`params.global.outdir`), by default ``out/``

- ``data``: contains the workflow output file (in h5ad format), plus symlinks to all the intermediate files.
- ``loom``: contains final loom files which can be imported inside SCope visualization tool for further visualization of the results.
- ``notebooks``: contains all the notebooks generated along the pipeline (e.g.: Quality control report)

    - See the example output report from the 1k PBMC data `here <http://htmlpreview.github.io/?https://github.com/vib-singlecell-nf/vsn-pipelines/blob/master/notebooks/10x_PBMC.merged_report.html>`_

- ``pipeline_reports``: nextflow dag, execution, timeline, and trace reports

If you would like to use the pipelines on a custom dataset, please see the `pipelines <./pipelines.html>`_ section below.
