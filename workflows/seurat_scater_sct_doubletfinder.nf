nextflow.preview.dsl=2

include {
	SEURAT__RNA_QC
} from '../src/seurat/processes/RNAQc.nf' params(params)

include {
	SEURAT__THRESHOLDFILTERING
} from '../src/seurat/processes/thresholdFiltering.nf' params(params)

include {
	SEURAT__SEURAT_TO_SCE
	SEURAT__SCE_TO_SEURAT_WITH_MERGE
} from '../src/seurat/processes/convertion.nf' params(params)

include {
	SINGLE_CELL_EXPERIMENT__NORMALIZATION
} from  '../src/singlecellexperiment/processes/normalization.nf' params(params)

include {
	SINGLE_CELL_EXPERIMENT__PCA_FILTERING
} from  '../src/singlecellexperiment/processes/pcaFiltering.nf' params(params)

include {
	SEURAT__SCTRANSFORM
} from  '../src/seurat/processes/scTransform.nf' params(params)

include {
	SEURAT__FIND_VARIABLE_FEATURES
} from '../src/seurat/processes/findVariableFeatures.nf' params(params)
include {
	SEURAT__MARKER_GENES
} from '../src/seurat/processes/markerGenes.nf' params(params)
include {
	SEURAT__DIMENSIONALITY_REDUCTION_PCA
} from '../src/seurat/processes/dimensionalityReduction.nf' params(params)
include {
	SEURAT__CLUSTERING
} from '../src/seurat/processes/clustering.nf' params(params)
include {
	SEURAT__ANNOTATION_GRAPHS
} from '../src/seurat/processes/annotationGraphs.nf' params(params)

include {
	DOUBLETFINDER__RUN
} from '../src/doubletfinder/processes/doubletFinder.nf' params(params)

workflow seurat_scater_sct_doubletfinder {
	take : seuratInput
	main :
		SEURAT__RNA_QC(seuratInput)
		SEURAT__THRESHOLDFILTERING(SEURAT__RNA_QC.out[0])
		SEURAT__SEURAT_TO_SCE(SEURAT__THRESHOLDFILTERING.out[0])
		SINGLE_CELL_EXPERIMENT__PCA_FILTERING(SEURAT__SEURAT_TO_SCE.out)
		SINGLE_CELL_EXPERIMENT__NORMALIZATION(SINGLE_CELL_EXPERIMENT__PCA_FILTERING.out[0])
		SEURAT__SCE_TO_SEURAT_WITH_MERGE(SINGLE_CELL_EXPERIMENT__NORMALIZATION.out,SEURAT__THRESHOLDFILTERING.out[0])
		SEURAT__SCTRANSFORM(SEURAT__SCE_TO_SEURAT_WITH_MERGE.out)
		SEURAT__FIND_VARIABLE_FEATURES(SEURAT__SCTRANSFORM.out,"SCT")
		SEURAT__MARKER_GENES(SEURAT__FIND_VARIABLE_FEATURES.out[0])
		SEURAT__DIMENSIONALITY_REDUCTION_PCA(SEURAT__MARKER_GENES.out[0],"SCT")
		SEURAT__CLUSTERING(SEURAT__DIMENSIONALITY_REDUCTION_PCA.out[0],"SCT")
		SEURAT__ANNOTATION_GRAPHS(SEURAT__CLUSTERING.out[0],"SCT")
		DOUBLETFINDER__RUN(SEURAT__CLUSTERING.out[0])

		if(params.R_rnaseq.outformat == "h5ad"){
			output = SEURAT__SEURAT_TO_ANNDATA(DOUBLETFINDER__RUN.out[0],"SCT")
		} else {
			output = DOUBLETFINDER__RUN.out[0]
		}

	emit:
		output
}
