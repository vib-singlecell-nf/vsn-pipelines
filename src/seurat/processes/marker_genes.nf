nextflow.enable.dsl=2

import java.nio.file.Paths
binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/seurat/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SEURAT__MARKER_GENES {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__MARKER_GENES.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.marker_genes)
        processParams = sampleParams.local
        """
        ${binDir}/deg/sc_marker_genes.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__MARKER_GENES.${processParams.off} \
            --seed ${params.global.seed} \
            ${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
            ${(processParams.containsKey('logfcTreshold')) ? '--logfc-threshold ' + processParams.logfcTreshold : ''} \
            ${(processParams.containsKey('minPct')) ? '--min-pct ' + processParams.minPct : ''} \
            ${(processParams.containsKey('onlyPos')) ? '--only-pos ' + processParams.onlyPos : ''}
        """
}
