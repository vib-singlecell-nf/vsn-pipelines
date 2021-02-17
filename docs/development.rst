Development Guide
=================

Create module
-------------

Tool-based modules are located in ``src/<tool-name>``, and each module has a specific structure for scripts and Nextflow processes (see `Repository structure`_ below).

Case study: Add `Harmony`
*************************

Harmony is a method published in `Nature Methods`_ that performs integration of single-cell data.

.. _`Nature Methods`: https://www.nature.com/articles/s41592-019-0619-0

Links:

- GitHub: https://github.com/immunogenomics/harmony
- Tutorial: https://github.com/immunogenomics/harmony/blob/master/vignettes/quickstart.Rmd


Steps:

#. Create a new issue on ``vsn-pipelines`` GitHub repository explaining which module you are going to add (e.g.: `Add Harmony batch correction method`).

#. `Fork the`_ ``vsn-pipelines`` repository to your own GitHub account (if you are an external collaborator).

    .. _`Fork the`: https://help.github.com/en/github/getting-started-with-github/fork-a-repo

#. From your local copy of ``vsn-pipelines`` GitHub repository, create a new branch called ``feature/[github-issue-id]-[description]``.

    In this case,

    - ``[github-issue-id] = 115``
    - ``[description] = add_harmony_batch_correction_method``

   It is highly recommended to start from the ``develop`` branch:

    .. code:: bash

        git checkout develop
        git fetch
        git pull
        git checkout -b feature/115-add_harmony_batch_correction_method

#. Use the `template repository`_ in the vib-singlecell-nf organisation to create the framework for the new module in ``src/<tool-name>``:

    .. code:: bash

        git clone --depth=1 https://github.com/vib-singlecell-nf/template.git src/harmony

    .. _`template repository`: https://github.com/vib-singlecell-nf/template

#. Now, you can start to edit file in the tool module that is now located in ``src/<tool-name>``.
   Optionally, you can delete the ``.git`` directory in the new module to avoid confusion in future local development:

    .. code:: bash

        rm -rf src/harmony/.git


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


#. Rename the ``nextflow.config`` file to create the ``harmony.config`` configuration file.

    * Each process's options should be in their own level. With a single process, you do not need one extra level.

    .. code:: groovy

        params {
            tools {
                harmony {
                    container = 'vibsinglecellnf/harmony:1.0'
                    report_ipynb = "${params.misc.test.enabled ? '../../..' : ''}/src/harmony/bin/reports/sc_harmony_report.ipynb"
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
        parser <- add_option(
        parser, 
        c("-s", "--seed"), 
        action = "store", 
        default = 617,
        help="Seed. [default %default]"
        )

        args <- parse_args(parser)

        cat("Parameters: \n")
        print(args)

        if(is.null(args$`vars-use`)) {
            stop("The parameter --vars-use has to be set.")
        }

        # Required by irlba::irlba (which harmony depends on) for reproducibility
        if(!is.null(args$seed)) {
        set.seed(args$seed)
        } else {
        warnings("No seed is set, this will likely give none reproducible results.")
        }

        input_ext <- tools::file_ext(args$`input-file`)

        if(input_ext == "h5ad") {
        # Current fix until https://github.com/satijalab/seurat/issues/2485 is fixed
        file <- hdf5r::h5file(filename = args$`input-file`, mode = 'r')
        if(!("X_pca" %in% names(x = file[["obsm"]]))) {
            stop("X_pca slot is not found in the AnnData (h5ad).")
        }
        obs <- file[['obs']][]
        pca_embeddings <- t(x = file[["obsm"]][["X_pca"]][,])
        row.names(x = pca_embeddings) <- obs$index
        colnames(x = pca_embeddings) <- paste0("PCA_", seq(from = 1, to = ncol(x = pca_embeddings)))
        metadata <- obs
        # seurat <- Seurat::ReadH5AD(file = args$`input-file`)
        # if(!("pca" %in% names(seurat@reductions)) || is.null(x = seurat@reductions$pca))
        #   stop("Expects a PCA embeddings data matrix but it does not exist.")
        # data <- seurat@reductions$pca
        # pca_embeddings <- data@cell.embeddings
        # metadata <- seurat@meta.data
        } else {
        stop(paste0("Unrecognized input file format: ", input_ext, "."))
        }

        print(paste0("PCA embeddings matrix has ", dim(x = data)[1], " rows, ", dim(x = data)[2], " columns."))

        if(sum(args$`vars-use` %in% colnames(x = metadata)) != length(x = args$`vars-use`)) {
            stop("Some argument value from the parameter(s) --vars-use are not found in the metadata.")
        }

        # Run Harmony
        # Expects PCA matrix (Cells as rows and PCs as columns.)
        harmony_embeddings <- harmony::HarmonyMatrix(
        data_mat = pca_embeddings
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



#. Create the Nextflow process that will run the Harmony R script defined in the previous step.

    .. code:: groovy

        nextflow.preview.dsl=2

        binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/harmony/bin/" : ""

        process SC__HARMONY__HARMONY_MATRIX {
            
            container params.tools.harmony.container
            publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
            clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=1:00:00 -A ${params.global.qsubaccount}"

            input:
                tuple val(sampleId), path(f)

            output:
                tuple val(sampleId), path("${sampleId}.SC__HARMONY__HARMONY_MATRIX.tsv")

            script:
                def sampleParams = params.parseConfig(sampleId, params.global, params.tools.harmony)
                processParams = sampleParams.local
                varsUseAsArguments = processParams.varsUse.collect({ '--vars-use' + ' ' + it }).join(' ')
                """
                ${binDir}run_harmony.R \
                    --seed ${params.global.seed} \
                    --input-file ${f} \
                    ${varsUseAsArguments} \
                    --output-prefix "${sampleId}.SC__HARMONY__HARMONY_MATRIX"
                """

        }


#. Create a Nextflow "subworkflow" that will call the Nextflow process defined in the previous step and perform some other tasks (dimensionality reduction, cluster identification, marker genes identification and report generation)

    This step is not required. However it this step is skipped, the code would still need to added into the main ``harmony`` workflow (`workflows/harmony.nf`, see the next step)

    .. code:: groovy

        nextflow.preview.dsl=2

        //////////////////////////////////////////////////////
        //  process imports:

        include {
            clean;
        } from '../../utils/processes/utils.nf' params(params)
        include {
            COMBINE_BY_PARAMS;
        } from "../../utils/workflows/utils.nf" params(params)
        include {
            PUBLISH as PUBLISH_BEC_OUTPUT;
            PUBLISH as PUBLISH_BEC_DIMRED_OUTPUT;
            PUBLISH as PUBLISH_FINAL_HARMONY_OUTPUT;
        } from "../../utils/workflows/utils.nf" params(params)

        include {
            SC__HARMONY__HARMONY_MATRIX;
        } from './../processes/runHarmony.nf' params(params)
        include {
            SC__H5AD_UPDATE_X_PCA;
        } from './../../utils/processes/h5adUpdate.nf' params(params)
        include {
            NEIGHBORHOOD_GRAPH;
        } from './../../scanpy/workflows/neighborhood_graph.nf' params(params)
        include {
            DIM_REDUCTION_TSNE_UMAP;
        } from './../../scanpy/workflows/dim_reduction.nf' params(params)
        include {
            SC__SCANPY__CLUSTERING_PARAMS;
        } from './../../scanpy/processes/cluster.nf' params(params)
        include {
            CLUSTER_IDENTIFICATION;
        } from './../../scanpy/workflows/cluster_identification.nf' params(params) // Don't only import a specific process (the function needs also to be imported)

        // reporting:
        include {
            GENERATE_DUAL_INPUT_REPORT
        } from './../../scanpy/workflows/create_report.nf' params(params)

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
                harmony_embeddings = SC__HARMONY__HARMONY_MATRIX( 
                    dimReductionData.map { 
                        it -> tuple(it[0], it[1])
                    } 
                )
                SC__H5AD_UPDATE_X_PCA( 
                    dimReductionData.map {
                        it -> tuple(it[0], it[1]) 
                    }.join(harmony_embeddings) 
                )

                PUBLISH_BEC_OUTPUT(
                    SC__H5AD_UPDATE_X_PCA.out,
                    "BEC_HARMONY.output",
                    "h5ad",
                    null,
                    false
                )

                NEIGHBORHOOD_GRAPH(
                    SC__H5AD_UPDATE_X_PCA.out.join(
                        dimReductionData.map { 
                            it -> tuple(it[0], it[2], *it[3..(it.size()-1)])
                        }
                    )
                )

                // Run dimensionality reduction
                DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )

                PUBLISH_BEC_DIMRED_OUTPUT(
                    DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
                    "BEC_HARMONY.dimred_output",
                    "h5ad",
                    null,
                    false
                )

                // Run clustering
                // Define the parameters for clustering
                def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.tools.scanpy.clustering) )
                CLUSTER_IDENTIFICATION(
                    normalizedTransformedData,
                    DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
                    "Post Batch Effect Correction (Harmony)"
                )

                marker_genes = CLUSTER_IDENTIFICATION.out.marker_genes.map {
                    it -> tuple(
                        it[0], // sampleId
                        it[1], // data
                        !clusteringParams.isParameterExplorationModeOn() ? null : it[2..(it.size()-1)], // Stash params
                    )
                }

                PUBLISH_FINAL_HARMONY_OUTPUT( 
                    marker_genes.map {
                        it -> tuple(it[0], it[1], it[2])
                    },
                    "BEC_HARMONY.final_output",
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
                    PUBLISH_FINAL_HARMONY_OUTPUT.out,
                    clusteringParams
                )
                harmony_report = GENERATE_DUAL_INPUT_REPORT(
                    becDualDataPrePost,
                    file(workflow.projectDir + params.tools.harmony.report_ipynb),
                    "SC_BEC_HARMONY_report",
                    clusteringParams.isParameterExplorationModeOn()
                )

            emit:
                data = CLUSTER_IDENTIFICATION.out.marker_genes
                cluster_report = CLUSTER_IDENTIFICATION.out.report
                harmony_report

        }

#. In the ``vsn-pipelines``, create a new main workflow called ``harmony.nf`` under ``workflows/``:

    .. code:: groovy

        nextflow.preview.dsl=2

        ////////////////////////////////////////////////////////
        //  Import sub-workflows/processes from the utils module:
        include {
            getBaseName
        } from '../src/utils/processes/files.nf'
        include {
            clean;
            SC__FILE_CONVERTER;
            SC__FILE_CONCATENATOR
        } from '../src/utils/processes/utils.nf' params(params)
        include {
            COMBINE_BY_PARAMS
        } from '../src/utils/workflows/utils.nf' params(params)
        include {
            SC__H5AD_TO_FILTERED_LOOM
        } from '../src/utils/processes/h5adToLoom.nf' params(params)
        include {
            FILE_CONVERTER
        } from '../src/utils/workflows/fileConverter.nf' params(params)
        include {
            UTILS__GENERATE_WORKFLOW_CONFIG_REPORT
        } from '../src/utils/processes/reports.nf' params(params)

        ////////////////////////////////////////////////////////
        //  Import sub-workflows/processes from the tool module:
        include {
            QC_FILTER
        } from '../src/scanpy/workflows/qc_filter.nf' params(params)
        include {
            NORMALIZE_TRANSFORM
        } from '../src/scanpy/workflows/normalize_transform.nf' params(params)
        include {
            HVG_SELECTION
        } from '../src/scanpy/workflows/hvg_selection.nf' params(params)
        include {
            NEIGHBORHOOD_GRAPH
        } from '../src/scanpy/workflows/neighborhood_graph.nf' params(params)
        include {
            DIM_REDUCTION_PCA
        } from '../src/scanpy/workflows/dim_reduction_pca.nf' params(params)
        include {
            DIM_REDUCTION_TSNE_UMAP
        } from '../src/scanpy/workflows/dim_reduction.nf' params(params)
        // cluster identification
        include {
            SC__SCANPY__CLUSTERING_PARAMS
        } from '../src/scanpy/processes/cluster.nf' params(params)
        include {
            CLUSTER_IDENTIFICATION
        } from '../src/scanpy/workflows/cluster_identification.nf' params(params)
        include {
            BEC_HARMONY
        } from '../src/harmony/workflows/bec_harmony.nf' params(params)
        // reporting:
        include {
            SC__SCANPY__MERGE_REPORTS
        } from '../src/scanpy/processes/reports.nf' params(params)
        include {
            SC__SCANPY__REPORT_TO_HTML
        } from '../src/scanpy/processes/reports.nf' params(params)


        workflow harmony {

            take:
                data

            main:
                out = data | \
                    SC__FILE_CONVERTER | \
                    FILTER_AND_ANNOTATE_AND_CLEAN

                if(params.tools.scanpy.containsKey("filter")) {
                    out = QC_FILTER( out ).filtered // Remove concat
                }
                if(params.utils.file_concatenator) {
                    out = SC__FILE_CONCATENATOR( 
                        out.map {
                            it -> it[1]
                        }.toSortedList( 
                            { a, b -> getBaseName(a, "SC") <=> getBaseName(b, "SC") }
                        ) 
                    )
                }
                if(params.tools.scanpy.containsKey("data_transformation") && params.tools.scanpy.containsKey("normalization")) {
                    out = NORMALIZE_TRANSFORM( out )
                }
                out = HVG_SELECTION( out )
                DIM_REDUCTION_PCA( out )
                NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )
                DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )

                // Perform the clustering step w/o batch effect correction (for comparison matter)
                clusterIdentificationPreBatchEffectCorrection = CLUSTER_IDENTIFICATION( 
                    NORMALIZE_TRANSFORM.out,
                    DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
                    "Pre Batch Effect Correction"
                )

                // Perform the batch effect correction
                BEC_HARMONY(
                    NORMALIZE_TRANSFORM.out,
                    // include only PCA since Harmony will correct this
                    DIM_REDUCTION_PCA.out,
                    clusterIdentificationPreBatchEffectCorrection.marker_genes
                )
                
                // Conversion
                // Convert h5ad to X (here we choose: loom format)
                if(params.utils?.file_concatenator) {
                    filteredloom = SC__H5AD_TO_FILTERED_LOOM( SC__FILE_CONCATENATOR.out )
                    scopeloom = FILE_CONVERTER(
                        BEC_HARMONY.out.data.groupTuple(),
                        'HARMONY.final_output',
                        'loom',
                        SC__FILE_CONCATENATOR.out
                    )
                } else {
                    filteredloom = SC__H5AD_TO_FILTERED_LOOM( SC__FILE_CONVERTER.out )
                    scopeloom = FILE_CONVERTER(
                        BEC_HARMONY.out.data.groupTuple(),
                        'HARMONY.final_output',
                        'loom',
                        SC__FILE_CONVERTER.out
                    )
                }
                
                project = CLUSTER_IDENTIFICATION.out.marker_genes.map { it -> it[0] }
                UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
                    file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
                )

                // Collect the reports:
                // Define the parameters for clustering
                def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.tools.scanpy.clustering) )
                // Pairing clustering reports with bec reports
                if(!clusteringParams.isParameterExplorationModeOn()) {
                    clusteringBECReports = BEC_HARMONY.out.cluster_report.map {
                        it -> tuple(it[0], it[1])
                    }.combine(
                        BEC_HARMONY.out.harmony_report.map {
                            it -> tuple(it[0], it[1])
                        },
                        by: 0
                    ).map {
                        it -> tuple(it[0], *it[1..it.size()-1], null)
                    }
                } else {
                    clusteringBECReports = COMBINE_BY_PARAMS(
                        BEC_HARMONY.out.cluster_report.map { 
                            it -> tuple(it[0], it[1], *it[2])
                        },
                        BEC_HARMONY.out.harmony_report,
                        clusteringParams
                    )
                }
                ipynbs = project.combine(
                    UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out
                ).join(
                    HVG_SELECTION.out.report.map {
                        it -> tuple(it[0], it[1])
                    }
                ).combine(
                    clusteringBECReports,
                    by: 0
                ).map {
                    it -> tuple(it[0], it[1..it.size()-2], it[it.size()-1])
                }

                // reporting:
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



#. Add a new Nextflow profile in the ``profiles`` section of the main ``nextflow.config`` of the ``vsn-pipelines`` repository:

    .. code:: groovy

        profiles {

            harmony {
                includeConfig 'src/scanpy/scanpy.config'
                includeConfig 'src/harmony/harmony.config'
            }
            ...
        }

#. Finally add a new entry in ``main.nf`` of the ``vsn-pipelines`` repository

    .. code:: groovy

        // run multi-sample with bbknn, output a scope loom file
        workflow harmony {

            include {
                harmony as HARMONY 
            } from './workflows/harmony' params(params)
            include {
                PUBLISH as PUBLISH_HARMONY 
            } from "./src/utils/workflows/utils" params(params)

            getDataChannel | HARMONY
            PUBLISH_HARMONY(
                HARMONY.out.scopeloom,
                "HARMONY",
                "loom",
                null,
                false
            )

        }

    You should now be able to configure (``nextflow config ...``) and run the ``harmony`` pipeline (``nextflow run ...``).

#. After confirming that your module is functional, you should create a pull request to merge your changes into the ``develop`` branch.

    - Make sure you have removed all references to ``TEMPLATE`` in your repository
    - Include some basic documentation for your module so people know what it does and how to use it.

   The pull request will be reviewed and accepted once it is confirmed to be working. Once the ``develop`` branch is merged into ``master``, the new tool will be part of the new release of VSN Pipelines!

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
            SC__CELLRANGER__MKFASTQ(file(params.tools.cellranger.mkfastq.csv), path(params.tools.cellranger.mkfastq.runFolder))
            SC__CELLRANGER__COUNT(file(params.tools.cellranger.count.transcriptome), SC__CELLRANGER__MKFASTQ.out.flatten())
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
        tools {
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
                container = 'docker://vib-singlecell-nf/scanpy:1.8.1'
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

