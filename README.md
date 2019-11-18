# SingleCellTxBenchmark

[![Nextflow](https://img.shields.io/badge/nextflow-19.10.0-brightgreen.svg)](https://www.nextflow.io/)

A repository of pipelines for single-cell data in Nextflow DSL2.

All the output generated by the pipelines will be located in the directory specified by `params.global.outdir` in `nextflow.config`.

# Dependencies

Make sure you have the following softwares installed,
- [Nextflow](https://www.nextflow.io/)
* A container system, either of:
    * [Docker](https://docs.docker.com/)
    * [Singularity](https://www.sylabs.io/singularity/)

# Quick start

To run a quick test of the single sample analysis pipeline, we can use the 1k PBMC datasets provided by 10x Genomics.
This will take only **~3min** to run.

1. The data first needs to be downloaded (instructions can be found 
    [here](data/README.md)).

2. Next, a config file needs to be generated.
In your working directory, run `nextflow config ...` with the appropriate profiles:
```bash
nextflow config aertslab/SingleCellTxBenchmark \
    -profile singularity,single_sample > single_sample.config
```
Now, edit `single_sample.config`.
Most of the default values are already set for the test dataset, but certain variables (e.g. container links) may need to be changed.
In particular, `params.global.tenx_folder` should point to the `outs/` folder in the 10x data, and
    `params.sc.file_converter` should be a path to the sample metadata file.

3. The pipeline can be run using the config file just generated (`-C ...`), and specifying the `single_sample` workflow as an entrypoint:
```bash
nextflow -C single_sample.config \
   run aertslab/SingleCellTxBenchmark \
      -entry single_sample
```

The pipelines will generate 3 types of results in the output directory (`params.global.outdir`) 
- `data`: contains all the intermediate files.
- `loom`: contains final loom files which can be imported inside SCope visualization tool for further insight of the results.
- `notebooks`: contains all the notebooks generated along the pipeline (e.g.: Quality control report)
    - See the example output report from the 1k PBMC data 
      [here](notebooks/10x_PBMC.merged_report.html).
- `pipeline_reports` (if `-profile report` was passed to `nextflow config ...`)

If you would like to use the pipelines on a custom dataset, please go to the `Pipelines` section (see below).

# Pipelines

## General workflow and strategy
### Running the pipeline directly from GitHub:
The intended usage for this pipeline is for the code to be run directly from GitHub.
This results in a separation of the Nextflow code and the results stored in the working directory.
For example:
```bash
nextflow run aertslab/SingleCellTxBenchmark \
    -profile singularity,single_sample \
    -entry single_sample
```
This picks up `aertslab/SingleCellTxBenchmark/main.nf` and runs the workflow defined by the `-entry` setting (here, `single_sample`), using the built-in configs, which are merged from each tool used (defined in the `single_sample` profile).
Specifying `nextflow run -latest ...` will download the latest commit prior to execution, or the `-r ...` option can be used to specify a specific commit or branch.
However, in nearly all cases it will be necessary to run the pipeline with a customized config file.

### Running the pipeline with a customized config file
The recommended method is to first run `nextflow config ...` to generate a complete config file (with the default parameters) in your working directory.
The tool-specific parameters, as well as Docker/Singularity profiles, are included when specifying the appropriate profiles to `nextflow config`.
Any of the parameters in this config file can then be edited and used to run the workflow of your choice.
For example, to run the `single_sample` workflow in a new working directory using the `singularity` profile:

1. Generate the config using the `single_sample` and `singularity` profiles:
```bash
mkdir single_sample_test && cd single_sample_test

nextflow config aertslab/SingleCellTxBenchmark \
    -profile singularity,single_sample > single_sample.config
```
2. Now run the workflow using the new config file (using `-C` to use **only** this file), specifying the proper workflow as the entry point:
```bash
nextflow -C single_sample.config \
   run aertslab/SingleCellTxBenchmark \
      -entry single_sample
```

## Single-sample Datasets
Pipelines to run a single sample or multiple samples separately.

### `single_sample`
![](https://github.com/aertslab/SingleCellTxBenchmark/workflows/single_sample/badge.svg)

The **single_sample** workflow will process 10x data,taking in 10x-structured data, and metadata file.
The standard analysis steps are run: filtering, normalization, log-transformation, HVG selection, dimensionality reduction, clustering, and loom file generation.
The output is a loom file with the results embedded.

### `single_sample_scenic`
![](https://github.com/aertslab/SingleCellTxBenchmark/workflows/single_sample_scenic/badge.svg)

Runs the `single_sample` workflow above, then runs the SCENIC workflow on the output, generating a comprehensive loom file with the combined results.
This could be very resource intensive, depending on the dataset.

### `scenic`
![](https://github.com/aertslab/SingleCellTxBenchmark/workflows/scenic/badge.svg)

Runs the SCENIC workflow alone, generating a loom file with only the SCENIC results.

### `nemesh`
Runs the Nemesh pipeline (Drop-seq) on a single sample or multiple samples separately.

Source: http://mccarrolllab.org/wp-content/uploads/2016/03/Drop-seqAlignmentCookbookv1.2Jan2016.pdf

### `scenic_multiruns`
Runs the SCENIC workflow multiple times (set by `params.sc.scenic.numRuns`), generating a loom file with the aggregated results from the multiple SCENIC runs.

## Multiple Datasets
Pipelines to aggregate multiple datasets together.

### `bbknn`
Runs the BBKNN pipeline (sample-specific filtering, merging of individual samples, normalization, log-transformation, HVG selection, PCA analysis, then the batch-effect correction steps: BBKNN, clustering, dimensionality reduction (UMAP only)).
The output is a loom file with the results embedded.

Source: https://github.com/Teichlab/bbknn/blob/master/examples/pancreas.ipynb

### `bbknn_scenic`
Runs the `bbknn` workflow above, then runs the SCENIC workflow on the output, generating a comprehensive loom file with the combined results.
This could be very resource intensive, depending on the dataset.


## Information on using 10xGenomics datasets

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

Setting the input directory appropriately will collect all the samples listed in the `filtered_feature_bc_matrix` directories listed above.
For example, in `params.global`, setting:
```
tenx_folder = "/home/data/cellranger/Sample*/outs/"
```
will recursively find all 10x samples in that directory.

# Repository structure

## Root
The repository root contains a `main.nf` and associated `nextflow.config`.
The root `main.nf` imports and calls sub-workflows defined in the modules.

## Modules
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
│   │   └── merge_SCENIC_motif_track_loom.py
│   ├── processes
│   │   ├── aucell.nf
│   │   ├── cistarget.nf
│   │   ├── grnboost2withoutDask.nf
│   │   └── mergeScenicLooms.nf
│   ├── main.nf
│   └── scenic.config
│
└── utils
    ├── bin
    │   ├── h5ad_to_loom.py
    │   ├── sc_file_annotator.py
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

## Workflows

Workflows (chains of nf processes) are defined in the module root folder (e.g. [src/Scanpy/bec_bbknn.nf](https://github.com/aertslab/SingleCellTxBenchmark/blob/module_refactor/src/scanpy/bec_bbknn.nf))
Workflows import multiple processes and define the workflow by name:
```groovy
include SC__CELLRANGER__MKFASTQ from './processes/mkfastq'  params(params)
include SC__CELLRANGER__COUNT   from './processes/count'    params(params)

workflow CELLRANGER {
    main:
        SC__CELLRANGER__MKFASTQ(file(params.sc.cellranger.mkfastq.csv), file(params.sc.cellranger.mkfastq.runFolder))
        SC__CELLRANGER__COUNT(file(params.sc.cellranger.count.transcriptome), SC__CELLRANGER__MKFASTQ.out.flatten())
    emit:
        SC__CELLRANGER__COUNT.out
}

```

### Workflow imports
Entire **sub-workflows** can also be imported in other workflows with one command (inheriting all of the process imports from the workflow definition):
```groovy
include CELLRANGER from '../cellranger/main.nf' params(params)
```

This leads to the ability to easily define **high-level workflows** in the master nf file: `aertslab/SingleCellTxBenchmark/main.nf`:
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


## Parameters structure
Parameters are stored in a separate config file per workflow, plus the main `nextflow.config`. 
These parameters are merged when starting the run using e.g.:
```groovy
includeConfig 'src/scenic/nextflow.config'
```

The parameter structure internally (post-merge) is:
```groovy
params {
    global {
        baseFilePath = "/opt/SingleCellTxBenchmark"
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
            container = 'docker://aertslab/sctx-scanpy:0.5.0'
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

# Development

## Module testing

Modules and processes can be tested independently, you can find an example in `src/utils/main.test.nf`.

The `SC__FILE_CONVERTER` process is tested against the `tiny` dataset available in `data/01.count`.

