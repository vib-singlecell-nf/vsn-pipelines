nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/directs/bin/" : ""

process SC__DIRECTS__SELECT_DEFAULT_CLUSTERING {

    container params.getToolParams("directs").container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
    label 'compute_resources__default'

    input:
        tuple \
            val(sampleId), \
            path(f), \
            val(stashedParams)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__DIRECTS__SELECT_DEFAULT_CLUSTERING.loom"), \
            val(stashedParams)

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("directs").select_default_clustering)
		processParams = sampleParams.local
        """
        ${binDir}select_default_clustering.py \
            ${f} \
            ${sampleId}.SC__DIRECTS__SELECT_DEFAULT_CLUSTERING.loom \
            ${(processParams.containsKey('cellEmbeddingsIndex')) ? '--cell-embeddings-index ' + processParams.cellEmbeddingsIndex : ''} \
            ${(processParams.containsKey('fromMinClusterSize')) ? '--from-min-cluster-size ' + processParams.fromMinClusterSize : ''} \
            ${(processParams.containsKey('toMinClusterSize')) ? '--to-min-cluster-size ' + processParams.toMinClusterSize : ''} \
            ${(processParams.containsKey('byMinClusterSize')) ? '--by-min-cluster-size ' + processParams.byMinClusterSize : ''} \
            ${(processParams.containsKey('fromMinSamples')) ? '--from-min-samples ' + processParams.fromMinSamples : ''} \
            ${(processParams.containsKey('toMinSamples')) ? '--to-min-samples ' + processParams.toMinSamples : ''} \
            ${(processParams.containsKey('byMinSamples')) ? '--by-min-samples ' + processParams.byMinSamples : ''}
        """
}

