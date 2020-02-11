Development
============

Repository structure
--------------------

Root
****

The repository root contains a ``main.nf`` and associated ``nextflow.config``.
The root ``main.nf`` imports and calls sub-workflows defined in the modules.

Modules
********
A "module" consists of a folder labeled with the tool name (Scanpy, SCENIC, utils, etc.), with subfolders for

* ``bin/`` (scripts passed into the container)
* ``processes/`` (where Nextflow processes are defined)

The root of the modules folder contains workflow files + associated configs (as many as there are workflows):

* ``main.nf`` + ``nextflow.config``
* ``single_sample.nf`` + ``scenic.config``
* ...

::

    src/
    ├── cellranger
    │   ├── main.nf
    │   ├── nextflow.config
    │   └── processes
    │       ├── count.nf
    │       └── mkfastq.nf
    │
    ├── channels
    │   └── tenx.nf
    │
    ├── scenic
    │   ├── bin
    │   │   ├── grnboost2_without_dask.py
    │   ├── processes
    │   │   ├── aucell.nf
    │   │   ├── cistarget.nf
    │   │   ├── grnboost2withoutDask.nf
    │   ├── main.nf
    │   └── scenic.config
    │
    └── utils
        ├── bin
        │   ├── h5ad_to_loom.py
        │   ├── sc_file_concatenator.py
        │   └── sc_file_converter.py
        ├── utils.config
        └── processes
            ├── files.nf
            ├── h5ad_to_loom.nf
            ├── utils_1.test.nf
            ├── utils_2.test.nf
            └── utils.nf

Workflows
*********

Workflows (chains of nf processes) are defined in the module root folder (e.g. `src/Scanpy/bec_bbknn.nf <https://github.com/vib-singlecell-nf/vsn-pipelines/blob/module_refactor/src/scanpy/bec_bbknn.nf>`_ )
Workflows import multiple processes and define the workflow by name:

.. code:: groovy

    include SC__CELLRANGER__MKFASTQ from './processes/mkfastq'  params(params)
    include SC__CELLRANGER__COUNT   from './processes/count'    params(params)

    workflow CELLRANGER {

        main:
            SC__CELLRANGER__MKFASTQ(file(params.sc.cellranger.mkfastq.csv), path(params.sc.cellranger.mkfastq.runFolder))
            SC__CELLRANGER__COUNT(file(params.sc.cellranger.count.transcriptome), SC__CELLRANGER__MKFASTQ.out.flatten())
        emit:
            SC__CELLRANGER__COUNT.out

    }


Workflow imports
****************

Entire **sub-workflows** can also be imported in other workflows with one command (inheriting all of the process imports from the workflow definition):

.. code:: groovy

    include CELLRANGER from '../cellranger/main.nf' params(params)

This leads to the ability to easily define **high-level workflows** in the master nf file: ``vib-singlecell-nf/vsn-pipelines/main.nf``:

.. code:: groovy

    include CELLRANGER from './src/cellranger/main.nf' params(params)
    include BEC_BBKNN from './src/scanpy/bec_bbknn.nf' params(params)
    include SCENIC from './src/scenic/main.nf' params(params)

    workflow {

        CELLRANGER()
        BEC_BBKNN( CELLRANGER.out )
        SCENIC( BEC_BBKNN.out )

    }

Parameters structure
********************

Parameters are stored in a separate config file per workflow, plus the main ``nextflow.config``.
These parameters are merged when starting the run using e.g.:

.. code:: groovy

    includeConfig 'src/scenic/nextflow.config'

The parameter structure internally (post-merge) is:

.. code:: groovy

    params {
        global {
            baseFilePath = "/opt/vib-singlecell-nf"
            project_name = "MCF7"
            ...
        }
        sc {
            utils {
                file_converter {
                    ...
                }
                file_annotator {
                    ...
                }
                file_concatenator {
                    ...
                }
            }
            scanpy {
                container = 'docker://vib-singlecell-nf/scanpy:0.5.0'
                filter {
                    ...
                }
                data_transformation {
                    ...
                }
                normalization {
                    ...
                }
                feature_selection {
                    ...
                }
                feature_scaling {
                    ...
                }
                dim_reduction {
                    pca {
                        dimReductionMethod = 'PCA'
                        ...
                    }
                    umap {
                        dimReductionMethod = 'UMAP'
                        ...
                    }
                }
                batch_effect_correct {
                    ...
                }
                clustering {
                    ...
                }
            }
        }
    }

Module testing
----------------

Modules and processes can be tested independently, you can find an example in ``src/utils/main.test.nf``.

The ``SC__FILE_CONVERTER`` process is tested against the ``tiny`` dataset available in ``data/01.count``.

