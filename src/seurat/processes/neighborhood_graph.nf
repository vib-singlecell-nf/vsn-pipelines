nextflow.enable.dsl=2

import java.nio.file.Paths
binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/seurat/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SEURAT__NEIGHBORHOOD_GRAPH {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)
    
    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__NEIGHBORHOOD_GRAPH.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.neighborhood_graph)
        processParams = sampleParams.local
        """
        ${binDir}/nn/sc_neighborhood_graph.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__NEIGHBORHOOD_GRAPH.${processParams.off} \
            ${(processParams.containsKey('nPcs')) ? '--n-pcs ' + processParams.nPcs : ''}
        """
}
