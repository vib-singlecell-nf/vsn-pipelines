nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/seurat/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SEURAT__NORMALIZATION {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)
    
    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__NORMALIZATION.${processParams.off}")
    
    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.normalization)
        processParams = sampleParams.local
        """
        ${binDir}/transform/sc_normalization.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__NORMALIZATION.${processParams.off} \
            ${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
            ${(processParams.containsKey('scaleFactor')) ? '--scale-factor ' + processParams.scaleFactor : ''}
        """
}

process SC__SEURAT__NORMALIZATION_SCT {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__NORMALIZATION_SCT.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sct)
        processParams = sampleParams.local
        """
        ${binDir}/transform/sc_normalization_SCT.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__NORMALIZATION_SCT.${processParams.off} \
            ${(processParams.containsKey('assayName')) ? '--new-assay ' + processParams.assayName : ''} \
            ${(processParams.containsKey('nCells')) ? '--n-cells ' + processParams.nCells : ''} \
            ${(processParams.containsKey('nVariableFeatures')) ? '--n-variable-features ' + processParams.nVariableFeatures : ''} \
            ${(processParams.containsKey('onlyVarFeatures')) ? '--only-var ' + processParams.onlyVarFeatures : ''} \
            --seed ${params.global.seed}
        """
}

process SC__SEURAT__SCALING {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__SCALING.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.normalization)
        processParams = sampleParams.local
        """
        ${binDir}/transform/sc_scaling.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__SCALING.${processParams.off} \
            ${(processParams.containsKey('onlyVar')) ? '--only-var ' + processParams.onlyVar : ''}
        """
}
