# Seurat Module

This module allows the analysis of scRNA-seq and CITE-seq data (with or without hashtags) from a count matrix using the R package Seurat

## Available process

The different steps are implemented in the "/processes" sub-directory on the shape of nextflow process.

## Workflows

Some small pipelines including multiple processes are pre-made in the "/workflows" directory.

It is possible to include processes and/or workflows in other nextflow scripts, taking into account the input format.

## Main script

Complete pipelines are present in the main.nf file.

## How to use

It's possible to use the module standalone or integrated in a bigger pipeline with several modules.

### Standalone uses

1. If the module is used with a container (Docker or Singularity) fill the configuration files :
	- Seurat/conf/container.config
	- Seurat/conf/standALoneRunContainer.config

2. Generate the configuration file using the right profiles among the followings :
	- technical profiles :
		- docker
		- singularity
		- standard
		- reports

	- pipeline profiles :
		- HTO
		- RNA
		- ADT
		- hashtags_citeseq
		- citeseq
		- hashtags_rnaseq

	- Run the nextflow command :
	 ```bash
	nextflow config nextflow.config -profile [chosen profiles],[chosen pipeline] > experimentName.config
	```
3. Fill the configuration file previously generated, with the following rules:
	- **parameterName = ""**, are mandatory
	- **parameterName = null**, are optional with default values specified in the tool manual
	- **parameterName = a value**, should not be changed for a default run

4. Enter the following command to run the desired pipeline :
	```bash
	nextflow -C experimentName.config run main.nf -entry selected pipeline
	```
