nextflow.enable.dsl=2

process PROCESS_SAMPLESHEET {
    publishDir "${params.global.outdir}/processed_samplesheet"

    input:
    path(sampleSheet)

    output:
    path("processedSampleSheet.csv"), emit: processedSampleSheet

    script:
    """
    sed -n '/Sample_Name/,\$p' ${sampleSheet} > processedSampleSheet.csv
    """

}