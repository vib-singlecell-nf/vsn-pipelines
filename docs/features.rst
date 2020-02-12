Features
=========

Change log fold change and FDR thresholds for markers stored in SCope loom
--------------------------------------------------------------------------

By default, log fold change and FDR thresholds are set to 0 and 0.05 respectively.
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

Select the optimal number of principal components
-------------------------------------------------

When generating the config using ``nextflow config`` (see above), add the ``pcacv`` profile.

Remarks:

- Make sure ``nComps`` config parameter (under ``dim_reduction`` > ``pca``) is not set.
- If ``nPcs`` is not set for t-SNE or UMAP config entries, then all the PCs from the PCA will be used in the computation.

Currently, only the Scanpy related pipelines have this feature implemented.

Cell-based metadata annotation
------------------------------

If you have (pre-computed) cell-based metadata and you'd like to add them as annotations, please read `cell-based metadata annotation <https://github.com/vib-singlecell-nf/vsn-pipelines/tree/develop/src/utils#cell-based-metadata-annotation>`_.

Sample-based metadata annotation
--------------------------------

If you have sample-based metadata and you'd like to annotate the cells with these annotations, please read `sample-based metadata annotation <https://github.com/vib-singlecell-nf/vsn-pipelines/tree/develop/src/utils#sample-based-metadata-annotation>`_.

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

The latest version only implements this feature for the following pipelines:

- ``single_sample``
- ``bbknn``

Since ``v0.9.0``, it is possible to explore several combinations of parameters. The current version (``v0.9.0``) of the VSN-Pipelines allows to explore the following parameters:

- ``params.sc.scanpy.clustering``

  - ``method`` ::

        clusteringMethods = ['louvain','leiden']

  - ``resolution`` ::

        resolutions = [0.4, 0.8]
