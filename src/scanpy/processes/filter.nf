nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SCANPY__COMPUTE_QC_STATS {

  	container params.tools.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

  	input:
        tuple val(sampleId), path(f)

	output:
        tuple val(sampleId), path("${sampleId}.SC__SCANPY__COMPUTE_QC_STATS.${processParams.off}")

	script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.scanpy.filter)
		processParams = sampleParams.local
        """
        ${binDir}/filter/sc_cell_gene_filtering.py \
            compute \
            $f \
            ${sampleId}.SC__SCANPY__COMPUTE_QC_STATS.${processParams.off} \
            ${processParams?.cellFilterStrategy ? '--cell-filter-strategy ' + processParams.cellFilterStrategy : ''} \
            ${processParams?.cellFilterMinNCounts ? '--min-n-counts ' + processParams.cellFilterMinNCounts : ''} \
            ${processParams?.cellFilterMaxNCounts ? '--max-n-counts ' + processParams.cellFilterMaxNCounts : ''} \
            ${processParams?.cellFilterMinNGenes ? '--min-n-genes ' + processParams.cellFilterMinNGenes : ''} \
            ${processParams?.cellFilterMaxNGenes ? '--max-n-genes ' + processParams.cellFilterMaxNGenes : ''} \
            ${processParams?.cellFilterMaxPercentMito ? '--max-percent-mito ' + processParams.cellFilterMaxPercentMito : ''} \
            ${processParams?.geneFilterMinNCells ? '--min-number-cells ' + processParams.geneFilterMinNCells : ''}
        """

}


process SC__SCANPY__GENE_FILTER {

    container params.tools.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

    input:
        tuple val(sampleId), path(f)

	output:
        tuple val(sampleId), path("${sampleId}.SC__SCANPY__GENE_FILTER.${processParams.off}")

	script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.scanpy.filter)
		processParams = sampleParams.local
        """
        ${binDir}/filter/sc_cell_gene_filtering.py \
            genefilter \
            $f \
            ${sampleId}.SC__SCANPY__GENE_FILTER.${processParams.off} \
            ${processParams?.geneFilterMinNCells ? '--min-number-cells ' + processParams.geneFilterMinNCells : ''}
        """

}


process SC__SCANPY__CELL_FILTER {

    container params.tools.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

    input:
        tuple val(sampleId), path(f)

	output:
        tuple val(sampleId), path("${sampleId}.SC__SCANPY__CELL_FILTER.${processParams.off}")
    
	script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.scanpy.filter)
		processParams = sampleParams.local
        """
        ${binDir}/filter/sc_cell_gene_filtering.py \
            cellfilter \
            $f \
            ${sampleId}.SC__SCANPY__CELL_FILTER.${processParams.off} \
            ${processParams?.cellFilterStrategy ? '--cell-filter-strategy ' + processParams.cellFilterStrategy : ''} \
            ${processParams?.cellFilterMinNCounts ? '--min-n-counts ' + processParams.cellFilterMinNCounts : ''} \
            ${processParams?.cellFilterMaxNCounts ? '--max-n-counts ' + processParams.cellFilterMaxNCounts : ''} \
            ${processParams?.cellFilterMinNGenes ? '--min-n-genes ' + processParams.cellFilterMinNGenes : ''} \
            ${processParams?.cellFilterMaxNGenes ? '--max-n-genes ' + processParams.cellFilterMaxNGenes : ''} \
            ${processParams?.cellFilterMaxPercentMito ? '--max-percent-mito ' + processParams.cellFilterMaxPercentMito : ''}
        """

}
