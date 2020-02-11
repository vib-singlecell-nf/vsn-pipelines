Getting Started
================

Dependencies
-------------
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
    Launching `vib-singlecell-nf/vsn-pipelines` [condescending_liskov] - revision: 92368248f3 [master]
    WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE

    [33/68d885] process > single_sample:SINGLE_SAMPLE:UTILS__GENERATE_WORKFLOW_CONFIG_REPORT                                          [100%] 1 of 1 ✔
    [a2/dcf990] process > single_sample:SINGLE_SAMPLE:QC_FILTER:SC__FILE_CONVERTER (1)                                                [100%] 1 of 1 ✔
    [9c/dff236] process > single_sample:SINGLE_SAMPLE:QC_FILTER:SC__SCANPY__COMPUTE_QC_STATS (1)                                      [100%] 1 of 1 ✔
    [65/e1bf9f] process > single_sample:SINGLE_SAMPLE:QC_FILTER:SC__SCANPY__GENE_FILTER (1)                                           [100%] 1 of 1 ✔
    [92/faae99] process > single_sample:SINGLE_SAMPLE:QC_FILTER:SC__SCANPY__CELL_FILTER (1)                                           [100%] 1 of 1 ✔
    [52/c39d90] process > single_sample:SINGLE_SAMPLE:QC_FILTER:GENERATE_DUAL_INPUT_REPORT:SC__SCANPY__GENERATE_DUAL_INPUT_REPORT (1) [100%] 1 of 1 ✔
    [d2/b38e10] process > single_sample:SINGLE_SAMPLE:QC_FILTER:GENERATE_DUAL_INPUT_REPORT:SC__SCANPY__REPORT_TO_HTML (1)             [100%] 1 of 1 ✔
    [87/96ef4d] process > single_sample:SINGLE_SAMPLE:NORMALIZE_TRANSFORM:SC__SCANPY__NORMALIZATION (1)                               [100%] 1 of 1 ✔
    [b2/493705] process > single_sample:SINGLE_SAMPLE:NORMALIZE_TRANSFORM:SC__SCANPY__DATA_TRANSFORMATION (1)                         [100%] 1 of 1 ✔
    [69/a2a237] process > single_sample:SINGLE_SAMPLE:HVG_SELECTION:SC__SCANPY__FEATURE_SELECTION (1)                                 [100%] 1 of 1 ✔
    [1d/0ec983] process > single_sample:SINGLE_SAMPLE:HVG_SELECTION:SC__SCANPY__FEATURE_SCALING (1)                                   [100%] 1 of 1 ✔
    [91/11965d] process > single_sample:SINGLE_SAMPLE:HVG_SELECTION:GENERATE_REPORT:SC__SCANPY__GENERATE_REPORT (1)                   [100%] 1 of 1 ✔
    [4e/620e9e] process > single_sample:SINGLE_SAMPLE:HVG_SELECTION:GENERATE_REPORT:SC__SCANPY__REPORT_TO_HTML (1)                    [100%] 1 of 1 ✔
    [fd/c6e8c5] process > single_sample:SINGLE_SAMPLE:DIM_REDUCTION:DIM_REDUCTION_PCA:SC__SCANPY__DIM_REDUCTION__PCA (1)              [100%] 1 of 1 ✔
    [32/548f80] process > single_sample:SINGLE_SAMPLE:DIM_REDUCTION:SC__SCANPY__DIM_REDUCTION__TSNE (1)                               [100%] 1 of 1 ✔
    [e0/9b68f3] process > single_sample:SINGLE_SAMPLE:DIM_REDUCTION:SC__SCANPY__DIM_REDUCTION__UMAP (1)                               [100%] 1 of 1 ✔
    [20/337908] process > single_sample:SINGLE_SAMPLE:DIM_REDUCTION:GENERATE_REPORT:SC__SCANPY__GENERATE_REPORT (1)                   [100%] 1 of 1 ✔
    [b9/dc2795] process > single_sample:SINGLE_SAMPLE:DIM_REDUCTION:GENERATE_REPORT:SC__SCANPY__REPORT_TO_HTML (1)                    [100%] 1 of 1 ✔
    [0b/42a0a3] process > single_sample:SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:SC__SCANPY__CLUSTERING (1)                               [100%] 1 of 1 ✔
    [3a/084e6f] process > single_sample:SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:GENERATE_REPORT:SC__SCANPY__GENERATE_REPORT (1)          [100%] 1 of 1 ✔
    [06/6ea130] process > single_sample:SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:GENERATE_REPORT:SC__SCANPY__REPORT_TO_HTML (1)           [100%] 1 of 1 ✔
    [84/ca1672] process > single_sample:SINGLE_SAMPLE:CLUSTER_IDENTIFICATION:SC__SCANPY__MARKER_GENES (1)                             [100%] 1 of 1 ✔
    [db/d66797] process > single_sample:SINGLE_SAMPLE:SC__H5AD_TO_FILTERED_LOOM (1)                                                   [100%] 1 of 1 ✔
    [46/be45d7] process > single_sample:SINGLE_SAMPLE:FILE_CONVERTER:SC__H5AD_TO_LOOM (1)                                             [100%] 1 of 1 ✔
    [78/3988ff] process > single_sample:SINGLE_SAMPLE:FILE_CONVERTER:COMPRESS_HDF5 (1)                                                [100%] 1 of 1 ✔
    [4d/bfb133] process > single_sample:SINGLE_SAMPLE:SC__PUBLISH_H5AD (1)                                                            [100%] 1 of 1 ✔
    [9c/b5f299] process > single_sample:SINGLE_SAMPLE:SC__SCANPY__MERGE_REPORTS (1)                                                   [100%] 1 of 1 ✔
    [00/b15be5] process > single_sample:SINGLE_SAMPLE:SC__SCANPY__REPORT_TO_HTML (1)                                                  [100%] 1 of 1 ✔
    Converting 1k_pbmc_v2_chemistry.SC__SCANPY__MARKER_GENES.h5ad to 1k_pbmc_v2_chemistry.SC__SCANPY__MARKER_GENES.loom (w/ additional compression)...
    Completed at: 22-Jan-2020 13:45:59
    Duration    : 2m 38s
    CPU hours   : 0.1
    Succeeded   : 28


The pipelines will generate 3 types of results in the output directory (`params.global.outdir`), by default `out/`

- ``data``: contains the workflow output file (in h5ad format), plus symlinks to all the intermediate files.
- ``loom``: contains final loom files which can be imported inside SCope visualization tool for further visualization of the results.
- ``notebooks``: contains all the notebooks generated along the pipeline (e.g.: Quality control report)

    - See the example output report from the 1k PBMC data `here <http://htmlpreview.github.io/?https://github.com/vib-singlecell-nf/vsn-pipelines/blob/master/notebooks/10x_PBMC.merged_report.html>`_

- ``pipeline_reports``: nextflow dag, execution, timeline, and trace reports

If you would like to use the pipelines on a custom dataset, please see the `pipelines <./pipelines.html>`_ section below.
