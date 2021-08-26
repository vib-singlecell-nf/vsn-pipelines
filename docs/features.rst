Advanced Features
=================

Two-pass strategy
-----------------

Typically, cell- and gene-level filtering is one of the first steps performed in the analysis pipelines.
This usually results in the pipeline being run in two passes.
In the **first pass**, the default filters are applied (which are probably not valid for new datasets), and a separate QC report is generated for each sample.
These QC reports can be inspected and the filters can be adjusted in the config file either for all samples (by editing the ``params.tools.scanpy.filter`` settings directly, or for individual samples by using the strategy described in multi-sample parameters.
Then, the **second pass** restarts the pipeline with the correct filtering parameters applied (use ``nextflow run ... -resume`` to skip already completed steps).

Other notes
^^^^^^^^^^^
In order to run a specific pipeline (e.g. ``single_sample``),
the pipeline name must be specified as a **profile** when running ``nextflow config ...`` (so that the default parameters are included),
and as the **entry** workflow when running the pipeline with ``nextflow run``.

One exception to this is that the ``-entry`` pipeline can be one that is a subset of the one present in the config file.
For example, in a pipeline with long running step that occurs after filtering (e.g. ``single_sample_scenic``),
it can be useful to generate the full config file (``nextflow config vib-singlecell-nf/vsn-pipelines -profile single_sample_scenic``),
then run a first pass for filtering using ``nextflow run vib-singlecell-nf/vsn-pipelines -entry single_sample``, and a second pass using the full pipeline ``-entry single_sample_scenic``).


Avoid re-running SCENIC and use pre-existing results
----------------------------------------------------
Often one would like to test different batch effect correction methods with SCENIC. Naively, one would run the following commands:

.. code:: bash

    nextflow config ~/vibsinglecellnf -profile tenx,bbknn,dm6,scenic,scenic_use_cistarget_motifs,singularity > bbknn.config
    nextflow -C bbknn.config run vib-singlecell-nf/vsn-pipelines -entry bbknn_scenic

and,

.. code:: bash

    nextflow config ~/vibsinglecellnf -profile tenx,bbknn,dm6,scenic,scenic_use_cistarget_motifs,singularity > bbknn.config
    nextflow -C harmony.config run vib-singlecell-nf/vsn-pipelines -entry harmony_scenic

The annoying bit here is that we run SCENIC twice. This is what we would like to avoid since the SCENIC results will be the same.
To avoid this one can run the following code for generating the `harmony_scenic.config`,

.. code:: bash

    nextflow config ~/vibsinglecellnf -profile tenx,harmony,scenic_append_only,singularity > harmony.config

This will add a different scenic entry in the config:

.. code:: bash

    params {
        tools {
            scenic {
                container = 'vibsinglecellnf/scenic:0.11.2'
                report_ipynb = '/src/scenic/bin/reports/scenic_report.ipynb'
                existingScenicLoom = ''
                sampleSuffixWithExtension = '' // Suffix after the sample name in the file path
                scenicoutdir = "${params.global.outdir}/scenic/"
                scenicScopeOutputLoom = 'SCENIC_SCope_output.loom'
            }
        }
    }

Make sure that the following entries are correctly set before running the pipeline,

- ``existingScenicLoom = ''``
- ``sampleSuffixWithExtension = '' // Suffix after the sample name in the file path``

Finally run the pipeline,

.. code:: bash

    nextflow -C harmony.config run vib-singlecell-nf/vsn-pipelines -entry harmony_scenic


Set the seed
------------
Some steps in the pipelines are non-deterministic. In order to have reproducible results, a seed is set by default to:

.. code:: groovy

    workflow.manifest.version.replaceAll("\\.","").toInteger()

The seed is a number derived from the version of the pipeline used at the time of the analysis run.
To override the seed (integer) you have edit the ``nextflow.config`` file with:

.. code:: groovy

    params {
        global {
            seed = [your-custom-seed]
        }
    }

This filter will only be applied on the final loom file of the VSN-Pipelines. All the intermediate files prior to the loom file will still contain all of them the markers.

Change log fold change (logFC) and false discovery rate (FDR) thresholds for the marker genes stored in the final SCope loom
----------------------------------------------------------------------------------------------------------------------------

By default, the logFC and FDR thresholds are set to 0 and 0.05 respectively.
If you want to change those thresholds applied on the markers genes, edit the ``nextflow.config`` with the following entries,

.. code:: groovy

    params {
        tools {
            scope {
                markers {
                    log_fc_threshold = 0.5
                    fdr_fc_threshold = 0.01
                }
            }
        }
    }

This filter will only be applied on the final loom file of the VSN-Pipelines. All the intermediate files prior to the loom file will still contain all of them the markers.

Automated selection of the optimal number of principal components
-----------------------------------------------------------------

When generating the config using ``nextflow config`` (see above), add the ``pcacv`` profile.

Remarks:

- Make sure ``nComps`` config parameter (under ``dim_reduction.pca``) is not set.
- If ``nPcs`` is not set for t-SNE or UMAP config entries, then all the PCs from the PCA will be used in the computation.

Currently, only the Scanpy related pipelines have this feature implemented.


.. _Cell annotation:

Cell-based metadata annotation
------------------------------

There are 2 ways of using this feature: either when running an end-to-end pipeline (e.g.: ``single_sample``, ``harmony``, ``bbknn``, ...) or on its own as a independent workflow.

Part of an end-to-end pipeline
******************************

The profile ``utils_cell_annotate`` should be added along with the other profiles when generating the main config using the ``nextflow config`` command.

For more detailed information about those parameters, please check the `cell_annotate parameter details <Parameters of cell_annotate_>`_ section below.

As an independent workflow
**************************

Please check the `cell_annotate`_ workflow.

.. _`cell_annotate`: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#nemesh

Parameters of cell_annotate
***************************

The ``utils_cell_annotate`` profile is adding the following part to the config:

.. code:: groovy

    params {
        tools {
            cell_annotate {
                off = 'h5ad'
                method = ''
                cellMetaDataFilePath = ''
                sampleSuffixWithExtension = ''
                indexColumnName = ''
                sampleColumnName = ''
                annotationColumnNames = ['']
            }
        }
    }

Two methods (``params.utils.cell_annotate.method``) are available:

- ``aio``
- ``obo``

If you have a single file containing the metadata information of all your samples, use ``aio`` method otherwise use ``obo``.

For both methods, here are the mandatory parameters to set:

- ``off`` should be set to ``h5ad``
- ``method`` choose either ``obo`` or ``aio``
- ``annotationColumnNames`` is an array of columns names from ``cellMetaDataFilePath`` containing different annotation metadata to add.

If ``aio`` used, the following additional parameters are required:

- ``cellMetaDataFilePath`` is a file path pointing to a single .tsv file (with header) with at least 2 columns: a column containing all the cell IDs and an annotation column.
- ``indexColumnName`` is the column name from ``cellMetaDataFilePath`` containing the cell IDs information. This column **can** have unique values; if it's not the case, it's important that the combination of the values from the ``indexColumnName`` and the ``sampleColumnName`` are unique. 
- ``sampleColumnName`` is the column name from ``cellMetaDataFilePath`` containing the sample ID/name information. Make sure that the values from this column match the samples IDs inferred from the data files. To know how those are inferred, please read the `Input Data Formats`_ section.

If ``obo`` is used, the following parameters are required:

- ``cellMetaDataFilePath``

  - In multi-sample mode, is a file path containing a glob pattern. The target file paths should each pointing to a .tsv file (with header) with at least 2 columns: a column containing all the cell IDs and an annotation column.
  - In single-sample mode, is a file path pointing to a single .tsv file (with header) with at least 2 columns: a column containing all the cell IDs and an annotation column.
  - **Note**: the file name(s) of ``cellMetaDataFilePath`` is/are required to contain the sample ID(s).

- ``sampleSuffixWithExtension`` is the suffix used to extract the sample ID from the file name(s) of ``cellMetaDataFilePath``. The suffix should be the part after the sample name in the file path.
- ``indexColumnName`` is the column name from ``cellMetaDataFilePath`` containing the cell IDs information. This column **must** have unique values.

.. _`Input Data Formats`: https://vsn-pipelines.readthedocs.io/en/develop/pipelines.html#input-data-formats


.. _Sample annotation:

Sample-based metadata annotation
--------------------------------

The profile ``utils_sample_annotate`` should be added when generating the main config using ``nextflow config``. This will add the following entry in the config:

.. code:: groovy

    params {
        tools {
            sample_annotate {
                iff = '10x_cellranger_mex'
                off = 'h5ad' 
                type = 'sample' 
                metadataFilePath = 'data/10x/1k_pbmc/metadata.tsv'
            }
        }
    }

Then, the following parameters should be updated to use the module feature:

- ``metadataFilePath`` is a .tsv file (with header) with at least 2 columns where the first column need to match the sample IDs. Any other columns will be added as annotation in the final loom i.e.: all the cells related to their sample will get annotated with their given annotations.

.. list-table:: Sample-based Metadata Table
    :widths: 40 40 20
    :header-rows: 1

    *   - id
        - chemistry
        - ...
    *   - 1k_pbmc_v2_chemistry
        - v2
        - ...
    *   - 1k_pbmc_v3_chemistry
        - v3
        - ...

Sample-annotating the samples using this system will allow any user to query all the annotation using the SCope portal. This is especially relevant when samples needs to be compared across specific annotations (check compare tab with SCope).

Cell-based metadata filtering
-----------------------------

There are 2 ways of using this feature: either when running an end-to-end pipeline (e.g.: ``single_sample``, ``harmony``, ``bbknn``, ...) or on its own as a independent workflow.

The ``utils_cell_filter`` profile is required when generating the config file. This profile will add the following part:

.. code:: groovy

    params {
        tools {
            cell_filter {
                off = 'h5ad'
                method = ''
                filters = [
                    [
                        id: '',
                        sampleColumnName: '',
                        filterColumnName: '',
                        valuesToKeepFromFilterColumn: ['']
                    ]
                ]
            }
        }
    }

Part of an end-to-end pipeline
******************************

For more detailed information about the parameters to set in ``params.utils.cell_filter``, please check the `cell_filter parameter details <Parameters of cell_filter>`_ section below.

As an independent workflow
**************************

Please check the `cell_filter`_ workflow or `cell_annotate_filter`_ workflow to perform cell-based annotation and cell-based filtering sequentially.

.. _`cell_filter`: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#cell_filter
.. _`cell_annotate_filter`: https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#cell_annotate_filter

Parameters of cell_filter
*************************

Two methods (``params.utils.cell_filter.method``) are available:

- ``internal``
- ``external``

If you have a single file containing the metadata information of all your samples, use ``external`` method otherwise use ``internal``.

For both methods, here are the mandatory parameters to set:

- ``off`` should be set to ``h5ad``
- ``method`` choose either ``internal`` or ``external``
- ``filters`` is a List of Maps where each Map is required to have the following parameters:

  - ``id`` is a short identifier for the filter
  - ``valuesToKeepFromFilterColumn`` is array of values from the ``filterColumnName`` that should be kept (other values will be filtered out).

If ``internal`` used, the following additional parameters are required:

- ``filters`` is a List of Maps where each Map is required to have the following parameters:

  - ``sampleColumnName`` is the column name containing the sample ID/name information. It should exist in the ``obs`` column attribute of the h5ad.
  - ``filterColumnName`` is the column name that will be used to filter out cells.  It should exist in the ``obs`` column attribute of the h5ad.

If ``external`` used, the following additional parameters are required:

- ``filters`` is a List of Maps where each Map is required to have the following parameters:

  - ``cellMetaDataFilePath`` is a file path pointing to a single .tsv file (with header) with at least 3 columns: a column containing all the cell IDs, another containing the sample ID/name information, and a column to use for the filtering.
  - ``indexColumnName`` is the column name from ``cellMetaDataFilePath`` containing the cell IDs information. This column **must** have unique values. 
  - `optional` ``sampleColumnName`` is the column name from ``cellMetaDataFilePath`` containing the sample ID/name information. Make sure that the values from this column match the samples IDs inferred from the data files. To know how those are inferred, please read the `Input Data Formats`_ section.
  - `optional` ``filterColumnName`` is the column name from ``cellMetaDataFilePath`` which be used to filter out cells.


.. _Multi-sample parameters:

Multi-sample parameters
------------------------

It's possible to define custom parameters for the different samples. It's as easy as defining a hashmap in groovy or a dictionary-like structure in Python.
You'll just have to repeat the following structure for the parameters which you want to enable the multi-sample feature for

.. code:: groovy

    params {
        tools {
            scanpy {
            container = 'vibsinglecellnf/scanpy:1.8.1'
            filter {
                report_ipynb = '/src/scanpy/bin/reports/sc_filter_qc_report.ipynb'
                // Here we enable the multi-sample feature for the cellFilterMinNgenes parameter
                cellFilterMinNGenes = [
                    '1k_pbmc_v2_chemistry': 600,
                    '1k_pbmc_v3_chemistry': 800
                ]
                // cellFilterMaxNGenes will be set to 4000 for all the samples
                cellFilterMaxNGenes = 4000
                // Here we again enable the multi-sample feature for the cellFilterMaxPercentMito parameter
                cellFilterMaxPercentMito = [
                    '1k_pbmc_v2_chemistry': 0.15,
                    '1k_pbmc_v3_chemistry': 0.05
                ]
                // geneFilterMinNCells will be set to 3 for all the samples
                geneFilterMinNCells = 3
                iff = '10x_mtx'
                off = 'h5ad'
                outdir = 'out'
            }
        }
    }

If you want to apply custom parameters for some specific samples and have a "general" parameter for the rest of the samples, you should use the 'default' key as follows:

.. code:: groovy

    params {
        tools {
            scanpy {
            container = 'vibsinglecellnf/scanpy:1.8.1'
            filter {
                report_ipynb = '/src/scanpy/bin/reports/sc_filter_qc_report.ipynb'
                // Here we enable the multi-sample feature for the cellFilterMinNgenes parameter
                cellFilterMinNGenes = [
                    '1k_pbmc_v2_chemistry': 600,
                    'default': 800
                ]
                [...]
            }
        }
    }

Using this config, the parameter ``params.tools.scanpy.cellFilterMinNGenes`` will be applied with a threshold value of ``600`` to ``1k_pbmc_v2_chemistry``.  The rest of the samples will use the value ``800`` to filter the cells having less than that number of genes.
This strategy can be applied to any other parameter of the config.


Parameter exploration
----------------------

Since ``v0.9.0``, it is possible to explore several combinations of parameters. The latest version of the VSN-Pipelines allows to explore the following parameters:

- ``params.tools.scanpy.clustering``

  - ``method`` ::

        methods = ['louvain','leiden']

  - ``resolution`` ::

        resolutions = [0.4, 0.8]

Select default clustering
*************************

In case the parameter exploration mode is used within the ``params.tools.scanpy.clustering`` parameter, it will generated a range of different clusterings. 
For non-expert, it's often difficult to know which clustering to pick. It's however possible to use the ``DIRECTS`` module in order to select a default clustering. In order, to use 
this automated clustering selection method, add the ``directs`` profile when generating the main config using ``nextflow config``. The config will get populated with:

.. code:: groovy

    directs {
        container = 'vibsinglecellnf/directs:0.1.0'
        labels {
            processExecutor = 'local'
        }
        select_default_clustering {
            fromMinClusterSize = 5
            toMinClusterSize = 100
            byMinClusterSize = 5
        }
    }

Currently, only the Scanpy related pipelines have this feature implemented.

Regress out variables
---------------------

By default, don't regress any variable out. To enable this features, the ``scanpy_regress_out`` profile should be added when generating the main config using ``nextflow config``. This will add the following entry in the config:

.. code:: groovy

    params {
        tools {
            scanpy {
                regress_out {
                    variablesToRegressOut = []
                    off = 'h5ad'
                }
            }
        }
    }

Add any variable in ``variablesToRegressOut`` to regress out: e.g.: 'n_counts', 'percent_mito'.

Highly Variable Genes Selection
-------------------------------

This step is a wrapper around the `Scanpy` ``scanpy.pp.highly_variable_genes`` function and regarding the parameters used it is following the documentation available at `scanpy-pp-highly-variable-genes`_.
By default, it will use the ``seurat`` flavor to select variable genes and will also keep the same default values for the 4 different thresholds (as the documentation): ``min_mean``, ``max_mean``, ``min_disp``, ``max_disp``.

.. _`scanpy-pp-highly-variable-genes`: https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.highly_variable_genes.html#scanpy-pp-highly-variable-genes.

.. code:: groovy

    params {
        tools {
            scanpy {
                feature_selection {
                    report_ipynb = "${params.misc.test.enabled ? '../../..' : ''}/src/scanpy/bin/reports/sc_select_variable_genes_report.ipynb"
                    flavor = 'seurat'
                    minMean = 0.0125
                    maxMean = 3
                    minDisp = 0.5
                    off = 'h5ad'
                }
            }
        }
    }


Other flavors are available as ``cell_ranger`` and ``seurat_v3``. In order to use the ``seurat_v3`` flavor, one parameter is required to be specified: ``nTopGenes`` in the config file as follows:

.. code:: groovy

    params {
        tools {
            scanpy {
                feature_selection {
                    report_ipynb = "${params.misc.test.enabled ? '../../..' : ''}/src/scanpy/bin/reports/sc_select_variable_genes_report.ipynb"
                    flavor = 'seurat_v3'
                    nTopGenes = 2000
                    off = 'h5ad'
                }
            }
        }
    }



Skip steps
----------

By default, the pipelines are run from raw data (unfiltered data, not normalized).

If you have already performed an independent steps with another it's possible to skip some steps from the pipelines. Currently, here are the steps that can be skipped:
- ``Scanpy`` filtering
- ``Scanpy`` normalization

Skip Scanpy filtering step
**************************

In order to skip the Scanpy filtering step, we need to add 3 new profiles when generating the config:

- ``min``
- ``scanpy_data_transformation``
- ``scanpy_normalization``

The following command, will create a Nextflow config which the pipeline will understand and will not run the Scanpy filtering step:

.. code:: groovy

    nextflow config \
       ~/vib-singlecell-nf/vsn-pipelines \
       -profile min,[data-profile],scanpy_data_transformation,scanpy_normalization,[...],singularity \
       > nextflow.config

- ``[data-profile]``: Can be one of the different possible data profiles e.g.: ``h5ad``
- ``[...]``: Can be other profiles like ``bbknn``, ``harmony``, ``pcacv``, ...


Quiet mode
----------

By default, VSN will output some additional messages to the terminal, such as the global seed, and the names and paths of the samples detected by the input channel.
These messages can be suppressed by using the ``--quiet`` flag when starting the nextflow process:

.. code:: bash

    nextflow -C example.config run vib-singlecell-nf/vsn-pipelines -entry single_sample --quiet

