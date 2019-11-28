nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/utils/bin/"
} else {
    binDir = ""
}

process SRP_TO_METADATA {

    publishDir "${params.outdir}/metadata", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=1 -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        file(sraProjectId)
    output:
        file "${sraProjectId}_metadata.tsv"
    script:
        """
        ${binDir}sra_to_metadata.py \
            ${sraProjectId} \
            --output "${sraProjectId}_metadata.tsv"
        """

}
