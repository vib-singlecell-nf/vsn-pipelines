nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/directs/bin/" : ""

process SC__DIRECTS__SELECT_DEFAULT_CLUSTERING {

    container params.sc.directs.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'

    input:
        tuple \
            val(sampleId), \
            path(f)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__DIRECTS__SELECT_DEFAULT_CLUSTERING.loom")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.directs.select_default_clustering)
		processParams = sampleParams.local
        """
        ${binDir}select_default_clustering.py \
            ${f} \
            ${sampleId}.SC__DIRECTS__SELECT_DEFAULT_CLUSTERING.loom \
            ${(processParams.containsKey('fromMinClusterSize')) ? '--from-min-cluster-size ' + processParams.fromMinClusterSize : ''} \
            ${(processParams.containsKey('toMinClusterSize')) ? '--to-min-cluster-size ' + processParams.toMinClusterSize : ''} \
            ${(processParams.containsKey('byMinClusterSize')) ? '--by-min-cluster-size ' + processParams.byMinClusterSize : ''}
        """
}

