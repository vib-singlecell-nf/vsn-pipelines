#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include { SEURAT__RNA_QC } from '../processes/RNAQc.nf' params(params)
include { SEURAT__THRESHOLDFILTERING } from '../processes/thresholdFiltering.nf' params(params)
include { SEURAT__SCTRANSFORM } from '../processes/scTransform.nf' params(params)
include { SEURAT__FIND_VARIABLE_FEATURES } from '../processes/findVariableFeatures.nf' params(params)
include { SEURAT__MARKER_GENES } from '../processes/markerGenes.nf' params(params)
include { SEURAT__DIMENSIONALITY_REDUCTION_PCA } from '../processes/dimensionalityReduction.nf' params(params)
include { SEURAT__CLUSTERING } from '../processes/clustering.nf' params(params)
include { SEURAT__ANNOTATION_GRAPHS } from '../processes/annotationGraphs.nf' params(params)
include {
			SEURAT__SEURAT_TO_SCE
			SEURAT__SCE_TO_SEURAT_WITH_MERGE
			SEURAT__SEURAT_TO_ANNDATA
		} from '../processes/convertion.nf' params(params)

workflow run_RNA {
	take: inputtuple
	main:
		SEURAT__RNA_QC(inputtuple)
		SEURAT__THRESHOLDFILTERING(SEURAT__RNA_QC.out[0])
		SEURAT__SCTRANSFORM(SEURAT__THRESHOLDFILTERING.out[0])
		SEURAT__FIND_VARIABLE_FEATURES(SEURAT__SCTRANSFORM.out,"SCT")
		SEURAT__MARKER_GENES(SEURAT__FIND_VARIABLE_FEATURES.out[0])
		SEURAT__DIMENSIONALITY_REDUCTION_PCA(SEURAT__MARKER_GENES.out[0],"SCT")
		SEURAT__CLUSTERING(SEURAT__DIMENSIONALITY_REDUCTION_PCA.out[0],"SCT")
		SEURAT__ANNOTATION_GRAPHS(SEURAT__CLUSTERING.out[0],"SCT")
		output = SEURAT__CLUSTERING.out[0]
		if(params.Seurat.RNA.outputFormat == "h5ad"){
			output = SEURAT__SEURAT_TO_ANNDATA(output,"SCT")
		}

	emit:
		output
}
