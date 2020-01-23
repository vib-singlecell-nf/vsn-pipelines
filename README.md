# VSN-Pipelines

[![Nextflow](https://img.shields.io/badge/nextflow-19.12.0-brightgreen.svg)](https://www.nextflow.io/)

A repository of pipelines for single-cell data in Nextflow DSL2.

This main repo contains multiple workflows for analyzing single cell transcriptomics data, and depends on a number of tools, which are organized into submodules within the [VIB-Singlecell-NF](https://github.com/vib-singlecell-nf) organization.
Currently available workflows include:
* **Cell Ranger**: processes 10x Chromium data to align reads to generate an expression counts matrix.
* **DropSeq**: processes Drop-seq data from read alignment to expression counts.
* **Single sample workflows**: perform a "best practices" scRNA-seq analysis. Multiple samples can be run in parallel, treating each sample separately.
* **Multi-sample workflows**: perform a "best practices" scRNA-seq analysis on a merged and batch-corrected group of samples. Available batch correction methods include:
  * **BBKNN**
  * **mnnCorrect**
* **GRN inference**:
  * The [pySCENIC](https://github.com/aertslab/pySCENIC) implementation of the [SCENIC](https://aertslab.org/#scenic) workflow is integrated here and can be run in conjunction with any of the above workflows.

The output of each of the main workflows is a [loom](http://loompy.org/)-format file, which is ready for import into the interactive single-cell web visualization tool [SCope](http://scope.aertslab.org/).
In addition, data is also output in h5ad format, and reports are generated for the major pipeline steps.

# Dependencies

Make sure you have the following software installed,
- [Nextflow](https://www.nextflow.io/)
* A container system, either of:
    * [Docker](https://docs.docker.com/)
    * [Singularity](https://www.sylabs.io/singularity/)

# Quick start

To run a quick test of the single sample analysis pipeline, we can use the 1k PBMC datasets provided by 10x Genomics.
This will take only **~3min** to run.

1. The data first needs to be downloaded (instructions can be found
    [here](data/README.md)).

2. Next, update to the latest pipeline version:
```bash
nextflow pull vib-singlecell-nf/vsn-pipelines
```

3. Next, generate a config file using the standard settings for the test data,
and the appropriate profiles (e.g., replace `singularity` with `docker` if necessary):
```bash
nextflow config vib-singlecell-nf/vsn-pipelines \
    -profile tenx,singularity,single_sample > single_sample.config
```

4. The test pipeline can now be run using the config file just generated, specifying the `single_sample` workflow as an entrypoint:
```bash
nextflow -C single_sample.config \
   run vib-singlecell-nf/vsn-pipelines \
      -entry single_sample
```

<details>
    <summary>Click here to see the expected output</summary>
    <pre>
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
    </pre>

</details>


The pipelines will generate 3 types of results in the output directory (`params.global.outdir`), by default `out/`
- `data`: contains the workflow output file (in h5ad format), plus symlinks to all the intermediate files.
- `loom`: contains final loom files which can be imported inside SCope visualization tool for further visualization of the results.
- `notebooks`: contains all the notebooks generated along the pipeline (e.g.: Quality control report)
    - See the example output report from the 1k PBMC data
      [here](notebooks/10x_PBMC.merged_report.html).
- `pipeline_reports`: nextflow dag, execution, timeline, and trace reports

If you would like to use the pipelines on a custom dataset, please see the `Pipelines` section below.


# Pipelines

## Generating a config file and running the pipeline

This pipeline can be configured and run on custom data with a few steps.
The recommended method is to first run `nextflow config ...` to generate a complete config file (with the default parameters) in your working directory.
The tool-specific parameters, as well as Docker/Singularity profiles, are included when specifying the appropriate profiles to `nextflow config`.

1. First, update to the latest pipeline version (this will update the nextflow cache of the repository, typically located in `~/.nextflow/assets/vib-singlecell-nf/`):
```bash
nextflow pull vib-singlecell-nf/vsn-pipelines
```

2. Next, a config file needs to be generated.
This step will merge parameters from multiple profiles together to create a master config which specifies **all** parameters used by the pipeline.
In this example, these are `tenx` for the input data, `singularity` to use the Singularity system (replace with `docker` if necessary), and `single_sample` to load the defaults for the single sample pipeline.
In your working directory, run `nextflow config ...` with the appropriate profiles:
```bash
nextflow config vib-singlecell-nf/vsn-pipelines \
    -profile tenx,singularity,single_sample > single_sample.config
```

3. Now, edits can be made to `single_sample.config`.
Generally, the default values are acceptable to use for a first pass, but certain variables (input directory, etc.) need to be changed.

In particular, the following parameters are frequently modified in practice:
* `params.global.project_name`: a project name which will be included in some of the output file names.
* `params.data.tenx.cellranger_outs_dir_path`, which should point to the `outs/` folder generated by CellRanger (if using 10x data).
  See [Information on using 10x Genomics datasets](#information-on-using-10x-genomics-datasets) for additional info.
* Filtering parameters (`params.sc.scanpy.filter`): filtering parameters, which will be applied to all samples, can be set here: min/max genes, mitochondrial read fraction, and min cells.
  See [Multi-sample parameters](#multi-sample-parameters) for additional info on how to specify sample-specific parameters.
* Louvain cluster resolution: `params.sc.scanpy.clustering.resolution`.
* For cell- and sample-level annotations, see [here](src/utils/README.md) for additional info.

4.  Run the workflow using the new config file (using `-C` is recommended to use **only** this file), specifying the proper workflow as the entry point:
```bash
nextflow -C single_sample.config \
   run vib-singlecell-nf/vsn-pipelines \
      -entry single_sample
```


### Two-pass strategy
Typically, cell- and gene-level filtering is one of the first steps performed in the analysis pipelines.
This usually results in the pipeline being run in two passes.
In the **first pass**, the default filters are applied (which are probably not valid for new datasets), and a separate QC report is generated for each sample.
These QC reports can be inspected and the filters can be adjusted in the config file
either for all samples (by editing the `params.sc.scanpy.filter` settings directly, or for individual samples by using the strategy described in [multi-sample parameters](#multi-sample-parameters).
Then, the **second pass** restarts the pipeline with the correct filtering parameters applied (use `nextflow run ... -resume` to skip already completed steps).

### Other notes
In order to run a specific pipeline (e.g. `single_sample`),
the pipeline name must be specified as a **profile** when running `nextflow config ...` (so that the default parameters are included),
and as the **entry** workflow when running the pipeline with `nextflow run`.

One exception to this is that the `-entry` pipeline can be one that is a subset of the one present in the config file.
For example, in a pipeline with long running step that occurs after filtering (e.g. `single_sample_scenic`),
it can be useful to generate the full config file (`nextflow config vib-singlecell-nf/vsn-pipelines -profile single_sample_scenic`),
then run a first pass for filtering using `nextflow run vib-singlecell-nf/vsn-pipelines -entry single_sample`, and a second pass using the full pipeline `-entry single_sample_scenic`).

## Single-sample Pipelines
Pipelines to run on a single sample or multiple samples separately and in parallel.

### `single_sample`
![](https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/single_sample/badge.svg)

The **single_sample** workflow will process 10x data,taking in 10x-structured data, and metadata file.
The standard analysis steps are run: filtering, normalization, log-transformation, HVG selection, dimensionality reduction, clustering, and loom file generation.
The output is a loom file with the results embedded.

<details>
    <summary>Click here to see the DAG of the workflow</summary>

![Single-sample Workflow](./assets/images/single_sample.svg)

</details>

### `single_sample_scenic`
![](https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/single_sample_scenic/badge.svg)

Runs the `single_sample` workflow above, then runs the SCENIC workflow on the output, generating a comprehensive loom file with the combined results.
This could be very resource intensive, depending on the dataset.

<details>
    <summary>Click here to see the DAG of the workflow</summary>

![Single-sample SCENIC Workflow](./assets/images/single_sample_scenic.svg)

</details>

### `scenic`
![](https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/scenic/badge.svg)

Runs the SCENIC workflow alone, generating a loom file with only the SCENIC results.

<details>
    <summary>Click here to see the DAG of the workflow</summary>

![SCENIC Workflow](./assets/images/scenic.svg)

</details>

### `scenic_multiruns`
Runs the SCENIC workflow multiple times (set by `params.sc.scenic.numRuns`), generating a loom file with the aggregated results from the multiple SCENIC runs.

<details>
    <summary>Click here to see the DAG of the workflow</summary>

![SCENIC Multi-runs Workflow](./assets/images/scenic_multiruns.svg)

</details>

### `cellranger`
Runs the cellranger workflow (`makefastq`, then `count`).
Input parameters are specified within the config file:
* `params.sc.cellranger.mkfastq.csv`: path to the CSV samplesheet
* `params.sc.cellranger.mkfastq.runFolder`: path of Illumina BCL run folder
* `params.sc.cellranger.count.transcriptome`: path to the Cell Ranger compatible transcriptome reference

### `nemesh`
Runs the Nemesh pipeline (Drop-seq) on a single sample or multiple samples separately.

[Source](http://mccarrolllab.org/wp-content/uploads/2016/03/Drop-seqAlignmentCookbookv1.2Jan2016.pdf)

## Multiple Datasets
Pipelines to aggregate multiple datasets together.

### `bbknn`
![](https://github.com/vib-singlecell-nf/vsn-pipelines/workflows/bbknn/badge.svg)

Runs the BBKNN pipeline (sample-specific filtering, merging of individual samples, normalization, log-transformation, HVG selection, PCA analysis, then the batch-effect correction steps: BBKNN, clustering, dimensionality reduction (UMAP only)).
The output is a loom file with the results embedded.

Source: https://github.com/Teichlab/bbknn/blob/master/examples/pancreas.ipynb

<details>
    <summary>Click here to see the DAG of the workflow</summary>

![BBKNN Workflow](./assets/images/bbknn.svg)

</details>

### `bbknn_scenic`
Runs the `bbknn` workflow above, then runs the SCENIC workflow on the output, generating a comprehensive loom file with the combined results.
This could be very resource intensive, depending on the dataset.

<details>
    <summary>Click here to see the DAG of the workflow</summary>

![BBKNN SCENIC Workflow](./assets/images/bbknn_scenic.svg)

</details>

## Information on using 10x Genomics datasets

Let's say the file structure of your data looks like this,

```
/home/data/
└── cellranger
    ├── Sample A
    │   └── outs
    │       ├── filtered_feature_bc_matrix
    │       └── ...
    └── Sample_B
        └── outs
            ├── filtered_feature_bc_matrix
            └── ...
```

Setting the input directory appropriately will collect all the samples listed in the `filtered_[feature|gene]_bc_matrix` directories listed above.
For example, in `params.data.tenx`, setting:
```
cellranger_outs_dir_path = "/home/data/cellranger/Sample*/outs/"
```
will recursively find all 10x samples in that directory.

# Docs

## Define input data

Depending on the type of data you run the pipeline with, one or more appropriate profiles should be set when running `nextflow config`.

### Cell Ranger (10xGenomics)
```
-profiles tenx
```
In the generated .config file, make sur the `cellranger_outs_dir_path` parameter is set with the paths to the Cell Ranger `outs` folders:
```
[...]
tenx {
    cellranger_outs_dir_path = "data/10x/1k_pbmc/1k_pbmc_*/outs/"
}
[...]
```
- The `cellranger_outs_dir_path` parameter accepts glob patterns and also comma separated paths.

### H5AD (Scanpy)
```
-profiles h5ad
```

In the generated .config file, make sur the `file_paths` parameter is set with the paths to the `.h5ad` files:
```
[...]
h5ad {
    file_paths = "data/1k_pbmc_v*_chemistry_SUFFIX.SC__FILE_CONVERTER.h5ad"
    suffix = "_SUFFIX.SC__FILE_CONVERTER.h5ad"
}
[...]
```
- The `suffix` parameter is used to infer the sample name from the file paths.
- The `file_paths` accepts glob patterns and also comma separated paths.
Make sure that `sc.file_converter.iff` is set to `h5ad`.

## Select the optimal number of principal components

When generating the config using `nextflow config` (see above), add the `pcacv` profile.

Remarks:
- Make sure `nComps` config parameter (under `dim_reduction` > `pca`) is not set.
- If `nPcs` is not set for t-SNE or UMAP config entries, then all the PCs from the PCA will be used in the computation.

Currently, only the Scanpy related pipelines have this feature implemented.

## Cell-based metadata annotation

If you have (pre-computed) cell-based metadata and you'd like to add them as annotations, please read [cell-based metadata annotation](https://github.com/vib-singlecell-nf/vsn-pipelines/tree/develop/src/utils#cell-based-metadata-annotation).

## Sample-based metadata annotation

If you have sample-based metadata and you'd like to annotate the cells with these annotations, please read [sample-based metadata annotation](https://github.com/vib-singlecell-nf/vsn-pipelines/tree/develop/src/utils#sample-based-metadata-annotation).

## Multi-sample parameters

It's possible to define custom parameters for the different samples. It's as easy as defining a hashmap in groovy or a dictionary-like structure in Python.
You'll just have to repeat the following structure for the parameters which you want to enable the multi-sample feature for:

```
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
```

# Development

## Repository structure

### Root
The repository root contains a `main.nf` and associated `nextflow.config`.
The root `main.nf` imports and calls sub-workflows defined in the modules.

### Modules
A "module" consists of a folder labeled with the tool name (Scanpy, SCENIC, utils, etc.), with subfolders for
* `bin/` (scripts passed into the container)
* `processes/` (where Nextflow processes are defined)
The root of the modules folder contains workflow files + associated configs (as many as there are workflows):
* `main.nf` + `nextflow.config`
* `single_sample.nf` + `scenic.config`
* ...

```
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
```

### Workflows

Workflows (chains of nf processes) are defined in the module root folder (e.g. [src/Scanpy/bec_bbknn.nf](https://github.com/vib-singlecell-nf/vsn-pipelines/blob/module_refactor/src/scanpy/bec_bbknn.nf))
Workflows import multiple processes and define the workflow by name:
```groovy
include SC__CELLRANGER__MKFASTQ from './processes/mkfastq'  params(params)
include SC__CELLRANGER__COUNT   from './processes/count'    params(params)

workflow CELLRANGER {

    main:
        SC__CELLRANGER__MKFASTQ(file(params.sc.cellranger.mkfastq.csv), path(params.sc.cellranger.mkfastq.runFolder))
        SC__CELLRANGER__COUNT(file(params.sc.cellranger.count.transcriptome), SC__CELLRANGER__MKFASTQ.out.flatten())
    emit:
        SC__CELLRANGER__COUNT.out

}

```

#### Workflow imports
Entire **sub-workflows** can also be imported in other workflows with one command (inheriting all of the process imports from the workflow definition):
```groovy
include CELLRANGER from '../cellranger/main.nf' params(params)
```

This leads to the ability to easily define **high-level workflows** in the master nf file: `vib-singlecell-nf/vsn-pipelines/main.nf`:
```groovy
include CELLRANGER from './src/cellranger/main.nf' params(params)
include BEC_BBKNN from './src/scanpy/bec_bbknn.nf' params(params)
include SCENIC from './src/scenic/main.nf' params(params)

workflow {

    CELLRANGER()
    BEC_BBKNN( CELLRANGER.out )
    SCENIC( BEC_BBKNN.out )

}
```


### Parameters structure
Parameters are stored in a separate config file per workflow, plus the main `nextflow.config`.
These parameters are merged when starting the run using e.g.:
```groovy
includeConfig 'src/scenic/nextflow.config'
```

The parameter structure internally (post-merge) is:
```groovy
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

```

## Module testing

Modules and processes can be tested independently, you can find an example in `src/utils/main.test.nf`.

The `SC__FILE_CONVERTER` process is tested against the `tiny` dataset available in `data/01.count`.

