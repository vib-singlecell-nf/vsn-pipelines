nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/utils/bin/"
} else {
    binDir = ""
}

process SRP_TO_METADATA {

    container params.sc.utils.container
    publishDir "${params.outdir}/metadata", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=1 -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        val(sraProjectId)
        val(sampleFilters)
    output:
        file "${sraProjectId}_metadata.tsv"
    script:
        sampleFiltersAsArguments = sampleFilters.collect({ '--sample-filter' + ' ' + it }).join(' ')
        """
        ${binDir}sra_to_metadata.py \
            ${sraProjectId} \
            ${sampleFiltersAsArguments} \
            --output "${sraProjectId}_metadata.tsv"
        """

}
