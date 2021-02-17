nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

def toolParams = params.tools.scenic

process CONVERT_MULTI_RUNS_FEATURES_TO_REGULONS {

    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/multi_runs_cistarget/", mode: 'link', overwrite: true
    // This process requires a large amount of memory especially for big datasets (force to use bigmem node)
    // This process is quite slow (could take more than 1h for big datasets, so keep 24h for now)
    label 'compute_resources__scenic_multiruns_motifs2regulons'

    input:
        tuple val(sampleId), path(multiRunsAggrMotifEnrichmentTable), path(multiRunsAggrRegulonsFolder)
        val type

    output:
        tuple val(sampleId), path("multi_runs_regulons_${type}.pkl.gz")

    script:
        """
        ${binDir}convert_multi_runs_features_to_regulons.py \
            $multiRunsAggrMotifEnrichmentTable \
            $multiRunsAggrRegulonsFolder \
            -o "multi_runs_regulons_${type}.pkl.gz"
        """

}
