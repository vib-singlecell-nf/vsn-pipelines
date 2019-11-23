nextflow.preview.dsl=2

if(!params.containsKey("test")) {
	binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
	binDir = ""
}

toolParams = params.sc.scenic
processParams = params.sc.scenic.motifs_to_regulons

process CONVERT_MULTI_RUNS_FEATURES_TO_REGULONS {

    // Process will be submitted as job if toolParams.labels.processExecutor = 'qsub' (default)
    label toolParams.labels.processExecutor
    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/multi_runs_cistarget/", mode: 'link', overwrite: true
    // This process requires a large amount of memory especially for big datasets (force to use bigmem node)
    // This process is quite slow (could take more than 1h for big datasets, so keep 24h for now)
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=${processParams.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    tuple val(sampleId), file(multiRunsAggrMotifEnrichmentTable), file(multiRunsAggrRegulonsFolder)
    val type

    output:
    tuple val(sampleId), file("multi_runs_regulons_${type}.pkl.gz")

    """
    ${binDir}convert_multi_runs_features_to_regulons.py \
        $multiRunsAggrMotifEnrichmentTable \
        $multiRunsAggrRegulonsFolder \
        -o "multi_runs_regulons_${type}.pkl.gz"
    """

}
