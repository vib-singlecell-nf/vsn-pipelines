#!/usr/bin/env nextflow
nextflow.preview.dsl=2


include {run_RNA } from'./workflows/RNA.nf' params(params)

workflow RNA {
	main:
	if(params.Seurat.seuratObjBuilder.inputFile == null){
		seuratInput = Channel.fromPath(params.Seurat.inputRdsFile)
						.map{ file -> tuple(params.sampleName,file)}
	} else {
		include {
					SEURAT__SEURAT_OBJECT_BUILDER
				} from './processes/seuratObjBuilder/seuratObjBuilder.nf' params(params)

		input = Channel.fromPath(params.Seurat.seuratObjBuilder.inputFile)
						.map{ file -> tuple(params.sampleName,file)}
		seuratInput = SEURAT__SEURAT_OBJECT_BUILDER(input)
	}
	run_RNA(seuratInput)
	emit:
	run_RNA.out
}
