nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/sratoolkit/bin/"
} else {
    binDir = ""
}

toolParams = params.sratoolkit

process DOWNLOAD_FASTQS_FROM_SRA_ACC_ID {

    // Process will be submitted as job if toolParams.labels.processExecutor = 'qsub' (default)
    label "${toolParams.labels ? toolParams.labels.processExecutor : "local"}"
    errorStrategy 'retry'
    maxRetries 3
    container toolParams.container
    publishDir "${params.global.outdir}/data/raw/fastqs", mode: 'symlink', overwrite: true
    clusterOptions "-l nodes=1:ppn=20 -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    maxForks toolParams.downloadFastqs.maxForks

    input:
        tuple val(sraId), val(sampleId)
    
    output:
        tuple val(sraId), file("${sraId}*.fastq.gz")
    
    script:
        if(sampleId == null || sampleId.length() < 1) {
            throw new Exception("DOWNLOAD_FASTQS_FROM_SRA_ACC_ID: Sample ID is empty.")
        }
        """
        prefetch -v -p 1 ${sraId}
        fasterq-dump -S -v -p -e ${toolParams.downloadFastqs.threads} -O . ${sraId}
        pigz -p ${toolParams.downloadFastqs.threads} *.fastq
        """

}
