nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/sratoolkit/bin/"
} else {
    binDir = ""
}

process DOWNLOAD_FASTQS_FROM_SRA_ACC_ID {

    errorStrategy 'retry'
    maxRetries 3
    container params.sratoolkit.container
    publishDir "${params.global.outdir}/data/raw/fastqs", mode: 'symlink', overwrite: true
    clusterOptions "-l nodes=1:ppn=20 -l walltime=24:00:00 -A ${params.qsubaccount}"
    maxForks params.sratoolkit.downloadFastqs.maxForks

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
        fasterq-dump -S -v -p -e ${params.sratoolkit.downloadFastqs.threads} -O . ${sraId}
        pigz -p ${params.sratoolkit.downloadFastqs.threads} *.fastq
        """

}
