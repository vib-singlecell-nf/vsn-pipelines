#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__SEURAT_TO_ANNDATA {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	container params.Seurat.container
  input:
	tuple val(sampleId), file(seuratobj)
	val(assay)
  output:
	tuple val(sampleId),file("${sampleId}.SEURAT__SEURAT_TO_ANNDATA_*.h5ad")
  script:
	"""
	Rscript ${scriptDir}/seuratToAnnData.R --seuratObj ${seuratobj} \
	--output "${sampleId}.SEURAT__SEURAT_TO_ANNDATA_${assay}.h5ad" \
	--assay "${assay}"
	"""
}

process SEURAT__SCE_TO_SEURAT_WITH_MERGE {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	container params.Seurat.container
  input:
	tuple val(sampleId), file(sceobj)
	tuple val(sampleId), file(seuratobj)
  output:
	tuple val(sampleId),file("${sampleId}.SEURAT__SCE_TO_SEURAT_WITH_MERGE.rds")
  script:
	"""
	Rscript ${scriptDir}/sceToSeurat.R --sceObj ${sceobj} \
	--seuratObj ${seuratobj} \
	--output "${sampleId}.SEURAT__SCE_TO_SEURAT_WITH_MERGE.rds"
	"""
}

process SEURAT__SCE_TO_SEURAT {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	container params.Seurat.container
  input:
	tuple val(sampleId), file(sceobj)
  output:
	tuple val(sampleId),file("${sampleId}.SEURAT__SCE_TO_SEURAT.rds")
  script:
	"""
	Rscript ${scriptDir}/sceToSeurat.R --sceObj ${sceobj}
	--output "${sampleId}.SEURAT__SCE_TO_SEURAT.rds"
	"""
}

process SEURAT__SEURAT_TO_SCE {
	//publishDir "${params.global.outDir}/${params.global.runName}", mode: 'symlink'
	container params.Seurat.container
  input:
	tuple val(sampleId), file(seuratobj)
  output:
	tuple val(sampleId),file("${sampleId}.SEURAT__SEURAT_TO_SCE.rds")
  script:
	"""
	Rscript ${scriptDir}/seuratToSce.R --seuratObj ${seuratobj} \
	--output "${sampleId}.SEURAT__SEURAT_TO_SCE.rds"
	"""
}
