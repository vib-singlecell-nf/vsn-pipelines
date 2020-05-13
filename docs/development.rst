Development
============

Create module
-------------

Case study: Add `Harmony`
*************************

Harmony is a method published in `Nature Methods`_ that performs integration of single-cell data.

.. _`Nature Methods`: https://www.nature.com/articles/s41592-019-0619-0

Links:

- GitHub: https://github.com/immunogenomics/harmony
- Tutorial: https://github.com/immunogenomics/harmony/blob/master/vignettes/quickstart.Rmd


Steps:

#. Ask the `VIB-SingleCell-NF` administrators to create a new repository (in this case: ``harmony``) or create one on your GitHub account that could be brought into the `VIB-SingleCell-NF` organization.

    When using your own repo, you MUST start from the `template repository`_ in the vib-singlecell-nf organisation. Click the green "Use this template" button and provide a name for your new repo. Make sure the "Include all branches" checkbox is checked.

    .. _`template repository`: https://github.com/vib-singlecell-nf/template

#. Create a new issue on ``vsn-pipelines`` GitHub repository explaining which module you are going to add (e.g.: `Add Harmony batch correction method`).


#. `Fork the`_ ``vsn-pipelines`` repository to your own GitHub account.

    .. _`Fork the`: https://help.github.com/en/github/getting-started-with-github/fork-a-repo

#. From your ``vsn-pipelines`` GitHub repository, create a new branch called ``feature/[github-issue-id]-[description]``.

    In this case,

    - ``[github-issue-id] = 115``
    - ``[description] = add_harmony_batch_correction_method``

    .. code:: bash

        git checkout -b feature/115-add_harmony_batch_correction_method

#. From within the ``src`` directory of the ``vsn-pipelines`` repo, run the ``add_new_submodule.sh`` script.

    .. code:: bash

        ./add_new_submodule.sh [git-repo-url] -d

    ``[git-repo-url]`` = https://github.com/vib-singlecell-nf/harmony.git (Git Repository URL from `VSN-SingleCell-NF` or from your GitHub account)
    ``-d`` tracks the develop branch of the new repository, which is where you should work until the module is working.

    If you are using VSCode and you don't see the new submodule appearing in ``SOURCE CONTROL PROVIDERS``, open any file from ``src/harmony`` (e.g.: LICENSE)


#. Create the Dockerfile recipe

    .. code:: dockerfile

        FROM dweemx/sctx-seurat:3.1.2

        RUN apt-get -y update && \
            apt-get install -y libcurl4-openssl-dev libxml2-dev zlib1g-dev libhdf5-dev && \
            apt-get install -y libssl-dev && \
            # png.h: No such file or directory
            apt-get install -y libpng-dev && \ 
            R -e "install.packages('optparse')" && \
            R -e "devtools::install_github(repo = 'aertslab/SCopeLoomR')" && \
            R -e "devtools::install_github('immunogenomics/harmony')" && \
            # Need to run ps
            apt-get -y install procps && \
            apt-get -y install libxml2 && \
            # Clean
            rm -rf /tmp/* && \
            apt-get autoremove -y && \
            apt-get autoclean -y && \
            rm -rf /var/cache/apt/* && \
            rm -rf /var/lib/apt/lists/* && \
            apt-get clean


#. Update the ``nextflow.config`` file to create the ``harmony.config`` configuration file.

    * Each process's options should be in their own level. With a single proccess, you do not need one extra level.

    .. code:: dockerfile

        params {
            sc {
                harmony {
                    container = 'vibsinglecellnf/harmony:1.0'
                    report_ipynb = "/src/harmony/bin/reports/sc_harmony_report.ipynb"
                    varsUse = ['batch']
                }
            }
        }

    The ``report_ipynb`` Jupyter Notebook is available here_.

    .. _here: https://github.com/vib-singlecell-nf/harmony/blob/master/bin/reports/sc_harmony_report.ipynb

#. Create the R script to run Harmony

    .. code:: r

        #!/usr/bin/env Rscript

        print("##################################################")
        print("# Harmony: Algorithm for single cell integration #")
        print("##################################################")

        # Loading dependencies scripts

        library("optparse")
            parser <- OptionParser(
            prog = "run_harmony.R",
            description = "Scalable integration of single cell RNAseq data for batch correction and meta analysis"
        )
        parser <- add_option(
            parser,
            c("-i", "--input-file"),
            action = "store",
            default = NULL,
            help = "Input file [default]"
        )
        parser <- add_option(
            parser,
            c("-a", "--vars-use"),
            action = "store",
            default = NULL,
            help = "If meta_data is dataframe, this defined which variable(s) to remove (character vector)."
        )
        parser <- add_option(
            parser,
            c("-p", "--do-pca"),
            action = "store",
            default = FALSE,
            help = "Whether to perform PCA on input matrix."
        )
        parser <- add_option(
            parser,
            c("-o", "--output-prefix"),
            action = "store",
            default = "foo",
            help="Prefix path to save output files. [default %default]"
        )

        args <- parse_args(parser)

        cat("Parameters: \n")
        print(args)

        if(is.null(args$`vars-use`)) {
            stop("The parameter --vars-use has to be set.")
        }

        input_ext <- tools::file_ext(args$`input-file`)

        if(input_ext == "h5ad") {
            seurat <- Seurat::ReadH5AD(file = args$`input-file`)
            if(!("pca" %in% names(seurat@reductions)) || is.null(x = seurat@reductions$pca))
                stop("Expects a PCA embeddings data matrix but it does not exist.")
            data <- seurat@reductions$pca
            metadata <- seurat@meta.data
        } else {
            stop(paste0("Unrecognized input file format: ", input_ext, "."))
        }

        print(paste0("PCA embeddings matrix has ", dim(x = data)[1], " rows, ", dim(x = data)[2], " columns."))

        if(sum(args$`vars-use` %in% colnames(x = metadata)) != length(x = args$`vars-use`)) {
            stop("Some argument value from the parameter(s) --vars-use are not found in the metadata.")
        }

        # Run Harmony
        harmony_embeddings <- harmony::HarmonyMatrix(data_mat = data@cell.embeddings
                                                    , meta_data = metadata
                                                    , vars_use = args$`vars-use`
                                                    , do_pca = args$`do-pca`
                                                    , verbose = FALSE
        )

        # Save the results

        ## PCA corrected embeddings

        write.table(
            x = harmony_embeddings,
            file = paste0(args$`output-prefix`, ".tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = TRUE,
            col.names = NA
        )


#. Create the Nextflow process that will run the Harmony R script defined in 7.

    .. code:: groovy

        nextflow.preview.dsl=2

        binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/harmony/bin/" : ""

        process SC__HARMONY__HARMONY_MATRIX {

            container params.sc.harmony.container
            publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
            clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=1:00:00 -A ${params.global.qsubaccount}"

            input:
                tuple val(sampleId), path(f)

            output:
                tuple val(sampleId), path("${sampleId}.SC__HARMONY__HARMONY_MATRIX.tsv")

            script:
                def sampleParams = params.parseConfig(sampleId, params.global, params.sc.harmony)
                processParams = sampleParams.local
                varsUseAsArguments = processParams.varsUse.collect({ '--vars-use' + ' ' + it }).join(' ')
                """
                ${binDir}run_harmony.R \
                    --input-file ${f} \
                    ${varsUseAsArguments} \
                    --output-prefix "${sampleId}.SC__HARMONY__HARMONY_MATRIX"
                """

        }

#. Create a Nextflow module that will call the Nextflow process defined in 8. and perform some other tasks (dimensionality reduction, cluster identification, marker genes identification and report generation)

    This step is not required. However it this step is skipped, the code would still need to added into the main ``harmony`` workflow (`workflows/harmony.nf`, see step 10)

    .. code:: groovy

        nextflow.preview.dsl=2

        //////////////////////////////////////////////////////
        //  process imports:

        include '../../utils/processes/utils.nf' params(params)
        include '../../utils/workflows/utils.nf' params(params)

        include SC__HARMONY__HARMONY_MATRIX from './../processes/runHarmony.nf' params(params)
        include SC__H5AD_UPDATE_X_PCA from './../../utils/processes/h5adUpdate.nf' params(params)
        include DIM_REDUCTION_TSNE_UMAP from './../../scanpy/workflows/dim_reduction.nf' params(params)
        include './../../scanpy/processes/cluster.nf' params(params)
        include './../../scanpy/workflows/cluster_identification.nf' params(params) // Don't only import a specific process (the function needs also to be imported)

        // reporting:
        include GENERATE_DUAL_INPUT_REPORT from './../../scanpy/workflows/create_report.nf' params(params)

        //////////////////////////////////////////////////////
        //  Define the workflow 

        workflow BEC_HARMONY {

            take:
                normalizedTransformedData
                dimReductionData
                // Expects (sampleId, anndata)
                clusterIdentificationPreBatchEffectCorrection

            main:
                // Run Harmony
                harmony_embeddings = SC__HARMONY__HARMONY_MATRIX( dimReductionData )
                SC__H5AD_UPDATE_X_PCA( 
                    dimReductionData.join(harmony_embeddings) 
                )

                // Run dimensionality reduction
                DIM_REDUCTION_TSNE_UMAP( SC__H5AD_UPDATE_X_PCA.out )

                // Run clustering
                // Define the parameters for clustering
                def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.sc.scanpy.clustering) )
                CLUSTER_IDENTIFICATION(
                    normalizedTransformedData,
                    DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
                    "Post Batch Effect Correction (Harmony)"
                )

                PUBLISH( 
                    CLUSTER_IDENTIFICATION.out.marker_genes.map { it -> tuple(it[0], it[1]) },
                    "BEC_HARMONY.output",
                    "h5ad",
                    null,
                    clusteringParams.isParameterExplorationModeOn()
                )

                // This will generate a dual report with results from
                // - Pre batch effect correction
                // - Post batch effect correction
                becDualDataPrePost = COMBINE_BY_PARAMS(
                    clusterIdentificationPreBatchEffectCorrection,
                    // Use PUBLISH output to avoid "input file name collision"
                    PUBLISH.out,
                    clusteringParams
                )
                harmony_report = GENERATE_DUAL_INPUT_REPORT(
                    becDualDataPrePost,
                    file(workflow.projectDir + params.sc.harmony.report_ipynb),
                    "SC_BEC_HARMONY_report",
                    clusteringParams.isParameterExplorationModeOn()
                )

            emit:
                data = CLUSTER_IDENTIFICATION.out.marker_genes
                cluster_report = CLUSTER_IDENTIFICATION.out.report
                harmony_report

        }

#. In the ``vsn-pipelines``, create a new main workflow called ``harmony.nf`` under ``workflows``

    .. code:: groovy

        nextflow.preview.dsl=2

        //////////////////////////////////////////////////////
        //  Import sub-workflows from the modules:

        include '../src/utils/processes/utils.nf' params(params.sc.file_concatenator + params.global + params)

        include QC_FILTER from '../src/scanpy/workflows/qc_filter.nf' params(params)
        include NORMALIZE_TRANSFORM from '../src/scanpy/workflows/normalize_transform.nf' params(params + params.global)
        include HVG_SELECTION from '../src/scanpy/workflows/hvg_selection.nf' params(params + params.global)
        include DIM_REDUCTION from '../src/scanpy/workflows/dim_reduction.nf' params(params + params.global)
        // CLUSTER_IDENTIFICATION
        include '../src/scanpy/processes/cluster.nf' params(params + params.global)
        include '../src/scanpy/workflows/cluster_identification.nf' params(params + params.global) // Don't only import a specific process (the function needs also to be imported)
        include BEC_HARMONY from '../src/harmony/workflows/bec_harmony.nf' params(params)

        include SC__H5AD_TO_FILTERED_LOOM from '../src/utils/processes/h5adToLoom.nf' params(params + params.global)
        include FILE_CONVERTER from '../src/utils/workflows/fileConverter.nf' params(params)

        // data channel to start from 10x data:
        include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

        // reporting:
        include UTILS__GENERATE_WORKFLOW_CONFIG_REPORT from '../src/utils/processes/reports.nf' params(params)
        include SC__SCANPY__MERGE_REPORTS from '../src/scanpy/processes/reports.nf' params(params + params.global)
        include SC__SCANPY__REPORT_TO_HTML from '../src/scanpy/processes/reports.nf' params(params + params.global)


        workflow harmony {

            take:
                data

            main:
                // Run the pipeline
                QC_FILTER( data ) // Remove concat 
                SC__FILE_CONCATENATOR( QC_FILTER.out.filtered.map{it -> it[1]}.collect() )
                NORMALIZE_TRANSFORM( SC__FILE_CONCATENATOR.out )
                HVG_SELECTION( NORMALIZE_TRANSFORM.out )
                DIM_REDUCTION( HVG_SELECTION.out.scaled )

                // Perform the clustering step w/o batch effect correction (for comparison matter)
                clusterIdentificationPreBatchEffectCorrection = CLUSTER_IDENTIFICATION( 
                    NORMALIZE_TRANSFORM.out,
                    DIM_REDUCTION.out.dimred_pca_tsne_umap,
                    "Pre Batch Effect Correction"
                )

                // Perform the batch effect correction
                BEC_HARMONY(
                    NORMALIZE_TRANSFORM.out,
                    // include only PCA since Harmony will correct this
                    DIM_REDUCTION.out.dimred_pca.map { it -> tuple(it[0], it[1]) },
                    clusterIdentificationPreBatchEffectCorrection.marker_genes
                )
                
                // Conversion
                // Convert h5ad to X (here we choose: loom format)
                filteredloom = SC__H5AD_TO_FILTERED_LOOM( SC__FILE_CONCATENATOR.out )
                scopeloom = FILE_CONVERTER(
                    BEC_HARMONY.out.data.groupTuple(),
                    'HARMONY.final_output'
                    'loom',
                    SC__FILE_CONCATENATOR.out
                )

                project = CLUSTER_IDENTIFICATION.out.marker_genes.map { it -> it[0] }
                UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
                    file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
                )
                // collect the reports:
                ipynbs = project.combine(
                    UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out
                ).join(
                    HVG_SELECTION.out.report
                ).join(
                    BEC_HARMONY.out.cluster_report
                ).combine(
                    BEC_HARMONY.out.harmony_report,
                    by: 0
                ).map {
                    tuple( it[0], it.drop(1) )
                }
                // reporting:
                def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.sc.scanpy.clustering) )
                SC__SCANPY__MERGE_REPORTS(
                    ipynbs,
                    "merged_report",
                    clusteringParams.isParameterExplorationModeOn()
                )
                SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)

            emit:
                filteredloom
                scopeloom

        }


#. Add a new Nextflow profile in ``nextflow.config`` of the ``vsn-pipelines`` repository

    .. code:: groovy

        workflow harmony {

            include harmony as HARMONY from './workflows/harmony' params(params)
            getDataChannel | HARMONY

        }

#. Finally add a new entry in main.nf of the ``vsn-pipelines`` repository

    .. code:: groovy

        harmony {
            includeConfig 'src/scanpy/scanpy.config'
            includeConfig 'src/harmony/harmony.config'
        }

    You should now be able to configure (``nextflow config``) and run the ``harmony`` pipeline (``nextflow run``).

#. After confirming that your module is functional, you should merge your changes in the tool repo into the ``master`` branch.

    - Make sure you have removed all references to ``TEMPLATE`` in your repository
    - Include some basic documentation for your module so people know what it does and how to use it.

#. Once merged into ``master`` you should update the submodule in the ``vsn-pipelines`` repo to point to the correct branch

    .. code:: bash

        git submodule set-branch --default src/harmony

#. Finally, add your new and updated files alongside the updated ``.gitmodules`` file and ``src/harmony`` files to a new commit and submit a pull request on the ``vsn-pipelines`` repo to have your new module integrated.

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
                        method = 'pca'
                        ...
                    }
                    umap {
                        method = 'tsne'
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

