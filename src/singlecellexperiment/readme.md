# Single Cell Experiment object Repository

This repository allows the analyze of scRNA-seq data from a count matrix using the R packages dropletUtils, scran, scater.

```
.
├── main.nf
├── modules
│   ├── bcRankMetrics
│   │   ├── bcRankMetrics.nf
│   │   ├── bcRankMetrics.R
│   │   └── readme.md
│   ├── cellQCMetrics
│   │   ├── cellQCMetrics.nf
│   │   ├── cellQCMetrics.R
|	|	└── readme.md
│   ├── filteringEmptyDroplets
│   │   ├── filteringEmptyDroplets.nf
│   │   └── filteringEmptyDroplets.R
│   ├── filteringThresholds
│   │   ├── filteringThresholds.nf
│   │   └── filteringThresholds.R
│   ├── normalization
│   │   ├── normalization.nf
│   │   └── normalization.R
│   ├── reportGenerator
│   │   ├── reportGenerator.nf
│   │   └── reportGenerator.Rmd
│   └── sceObjBuilder
│       ├── sceBuilder.nf
│       ├── sceBuilder.R
|		└── readme.md
├── nextflow.config
├── readme.md
├── runtest.sh
└── workflows
    └── filtering.nf

```

## Available process

All analyses implemented are in the "/modules" sub-directory on the shape of nextflow process. For more information, check the readme associated to each module.

## Workflows

Some small pipelines including multiple processes are pre-made in the "/workflows" directory.

It is possible to include processes and/or workflows in other nextflow scripts, taking into account the input format.

## Main script

A pre-build pipeline is present in the main.nf file.
