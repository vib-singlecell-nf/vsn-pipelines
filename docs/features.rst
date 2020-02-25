Features
=========

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
        sc {
            scenic {
                container = 'vibsinglecellnf/scenic:0.9.19'
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
Some steps in the pipelines are nondeterministic. To be able that the results are reproducible in time, by default a seed is set to:

.. code:: groovy

    workflow.manifest.version.replaceAll("\\.","").toInteger()

The seed is a number derived from the the version of the pipeline used at the time of the analysis run.
To override the seed (integer) you have edit the nextflow.config file with:

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
        sc {
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

- Make sure ``nComps`` config parameter (under ``dim_reduction`` > ``pca``) is not set.
- If ``nPcs`` is not set for t-SNE or UMAP config entries, then all the PCs from the PCA will be used in the computation.

Currently, only the Scanpy related pipelines have this feature implemented.

Cell-based metadata annotation
------------------------------

The profile ``utils_cell_annotate`` should be added when generating the main config using ``nextflow config``. This will add the following entry in the config:

.. code:: groovy

    params {
        sc {
            cell_annotate {
                iff = '10x_cellranger_mex'
                off = 'h5ad'
                cellMetaDataFilePath = ''
                indexColumnName = ''
                sampleColumnName = ''
                annotationColumnNames = ['']
            }
        }
    }

Then, the following parameters should be updated to use the module feature:

- ``cellMetaDataFilePath`` is a TSV file (with header) with at least 2 columns: a column containing all the cell IDs and an annotation column.
- ``indexColumnName`` is the column name from ``cellMetaDataFilePath`` containing the cell IDs information.
- ``sampleColumnName`` is the column name from ``cellMetaDataFilePath`` containing the sample ID/name information.
- ``annotationColumnNames`` is an array of columns names from ``cellMetaDataFilePath`` containing different annotation metadata to add.

Sample-based metadata annotation
--------------------------------

The profile ``utils_sample_annotate`` should be added when generating the main config using ``nextflow config``. This will add the following entry in the config:

.. code:: groovy

    params {
        sc {
            sample_annotate {
                iff = '10x_cellranger_mex'
                off = 'h5ad' 
                type = 'sample' 
                metaDataFilePath = 'data/10x/1k_pbmc/metadata.tsv'
            }
        }
    }

Then, the following parameters should be updated to use the module feature:

- ``metaDataFilePath`` is a TSV file (with header) with at least 2 columns where the first column need to match the sample IDs. Any other columns will be added as annotation in the final loom i.e.: all the cells related to their sample will get annotated with their given annotations.

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

Multi-sample parameters
------------------------

It's possible to define custom parameters for the different samples. It's as easy as defining a hashmap in groovy or a dictionary-like structure in Python.
You'll just have to repeat the following structure for the parameters which you want to enable the multi-sample feature for

.. code:: groovy

    params {
        sc {
            scanpy {
            container = 'vibsinglecellnf/scanpy:0.5.0'
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

Parameter exploration
----------------------

Since ``v0.9.0``, it is possible to explore several combinations of parameters. The latest version of the VSN-Pipelines allows to explore the following parameters:

- ``params.sc.scanpy.clustering``

  - ``method`` ::

        methods = ['louvain','leiden']

  - ``resolution`` ::

        resolutions = [0.4, 0.8]

Regress out variables
---------------------

By default, don't regress any variable out. To enable this features, the ``scanpy_regress_out`` profile should be added when generating the main config using ``nextflow config``. This will add the following entry in the config:

.. code:: groovy

    params {
        sc {
            scanpy {
                regress_out {
                    variablesToRegressOut = []
                    off = 'h5ad'
                }
            }
        }
    }

Add any variable in ``variablesToRegressOut`` to regress out: e.g.: 'n_counts', 'percent_mito'.
