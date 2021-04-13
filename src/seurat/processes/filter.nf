nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/seurat/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SEURAT__FEATURE_FILTER {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__FEATURE_FILTER.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.filter)
        processParams = sampleParams.local
        """
        ${binDir}/filter/sc_cell_gene_filtering.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__FEATURE_FILTER.${processParams.off} \
            --type feature \
            ${processParams.containsKey('featureFilterMinNCells') ? '--min-number-cells ' + processParams.featureFilterMinNCells : ''}
        """

}

process SC__SEURAT__CELL_FILTER {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__CELL_FILTER.${processParams.off}")
    
    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.filter)
        processParams = sampleParams.local
        """
        ${binDir}/filter/sc_cell_gene_filtering.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__CELL_FILTER.${processParams.off} \
            --type cell \
            ${processParams.containsKey('cellFilterMinNCounts') ? '--min-n-counts ' + processParams.cellFilterMinNCounts : ''} \
            ${processParams.containsKey('cellFilterMaxNCounts') ? '--max-n-counts ' + processParams.cellFilterMaxNCounts : ''} \
            ${processParams.containsKey('cellFilterMinNFeatures') ? '--min-n-features ' + processParams.cellFilterMinNFeatures : ''} \
            ${processParams.containsKey('cellFilterMaxNFeatures') ? '--max-n-features ' + processParams.cellFilterMaxNFeatures : ''} \
            ${processParams.containsKey('cellFilterMaxPercentMito') ? '--max-percent-mito ' + processParams.cellFilterMaxPercentMito : ''} \
            ${processParams.containsKey('mitoPrefix') ? '--mito-prefix ' + processParams.mitoPrefix : ''}
        """

}
