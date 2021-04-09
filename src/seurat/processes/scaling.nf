nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/seurat/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SEURAT__SCALING {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)
    
    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__SCALING.${processParams.off}")
    
    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.scaling)
        processParams = sampleParams.local
        """
        ${binDir}/transform/sc_scaling.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__SCALING.${processParams.off}" \
        """
}