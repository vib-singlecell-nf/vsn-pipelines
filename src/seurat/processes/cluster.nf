nextflow.enable.dsl=2

import java.nio.file.Paths
binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/seurat/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SEURAT__CLUSTERING {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__CLUSTERING.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.clustering)
        processParams = sampleParams.local
        """
        ${binDir}/clustering/sc_clustering.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__CLUSTERING.${processParams.off} \
            --seed ${params.global.seed} \
            ${(processParams.containsKey('algorithm')) ? '--algorithm ' + processParams.algorithm : ''} \
            ${(processParams.containsKey('resolution')) ? '--resolution ' + processParams.resolution : ''}
        """
}
