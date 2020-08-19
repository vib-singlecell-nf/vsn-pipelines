nextflow.preview.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SCANPY__COMPUTE_QC_STATS {

  	container params.sc.scanpy.container
    label 'compute_resources__mem'

  	input:
        tuple val(sampleId), path(f)

	output:
        tuple val(sampleId), path("${sampleId}.SC__SCANPY__COMPUTE_QC_STATS.${processParams.off}")

	script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.filter)
		processParams = sampleParams.local
        """
        ${binDir}/filter/sc_cell_gene_filtering.py \
            compute \
            $f \
            ${sampleId}.SC__SCANPY__COMPUTE_QC_STATS.${processParams.off} \
            ${(processParams.containsKey('cellFilterMinNCounts')) ? '--min-n-counts ' + processParams.cellFilterMinNCounts : ''} \
            ${(processParams.containsKey('cellFilterMaxNCounts')) ? '--max-n-counts ' + processParams.cellFilterMaxNCounts : ''} \
            ${(processParams.containsKey('cellFilterMinNGenes')) ? '--min-n-genes ' + processParams.cellFilterMinNGenes : ''} \
            ${(processParams.containsKey('cellFilterMaxNGenes')) ? '--max-n-genes ' + processParams.cellFilterMaxNGenes : ''} \
            ${(processParams.containsKey('cellFilterMaxPercentMito')) ? '--max-percent-mito ' + processParams.cellFilterMaxPercentMito : ''} \
            ${(processParams.containsKey('geneFilterMinNCells')) ? '--min-number-cells ' + processParams.geneFilterMinNCells : ''}
        """

}


process SC__SCANPY__GENE_FILTER {

    container params.sc.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

    input:
        tuple val(sampleId), path(f)

	output:
        tuple val(sampleId), path("${sampleId}.SC__SCANPY__GENE_FILTER.${processParams.off}")

	script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.filter)
		processParams = sampleParams.local
        """
        ${binDir}/filter/sc_cell_gene_filtering.py \
            genefilter \
            $f \
            ${sampleId}.SC__SCANPY__GENE_FILTER.${processParams.off} \
            ${(processParams.containsKey('geneFilterMinNCells')) ? '--min-number-cells ' + processParams.geneFilterMinNCells : ''}
        """

}


process SC__SCANPY__CELL_FILTER {

    container params.sc.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

    input:
        tuple val(sampleId), path(f)

	output:
        tuple val(sampleId), path("${sampleId}.SC__SCANPY__CELL_FILTER.${processParams.off}")
    
	script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.filter)
		processParams = sampleParams.local
        """
        ${binDir}/filter/sc_cell_gene_filtering.py \
            cellfilter \
            $f \
            ${sampleId}.SC__SCANPY__CELL_FILTER.${processParams.off} \
            ${(processParams.containsKey('cellFilterMinNCounts')) ? '--min-n-counts ' + processParams.cellFilterMinNCounts : ''} \
            ${(processParams.containsKey('cellFilterMaxNCounts')) ? '--max-n-counts ' + processParams.cellFilterMaxNCounts : ''} \
            ${(processParams.containsKey('cellFilterMinNGenes')) ? '--min-n-genes ' + processParams.cellFilterMinNGenes : ''} \
            ${(processParams.containsKey('cellFilterMaxNGenes')) ? '--max-n-genes ' + processParams.cellFilterMaxNGenes : ''} \
            ${(processParams.containsKey('cellFilterMaxPercentMito')) ? '--max-percent-mito ' + processParams.cellFilterMaxPercentMito : ''}
        """

}
