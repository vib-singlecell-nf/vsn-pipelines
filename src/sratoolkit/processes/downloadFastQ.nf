nextflow.enable.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/sratoolkit/bin/"
} else {
    binDir = ""
}

toolParams = params.tools.sratoolkit

process DOWNLOAD_FASTQS_FROM_SRA_ACC_ID {

    container toolParams.container
    publishDir "${params.global.outdir}/data/raw/fastqs", mode: 'symlink', overwrite: true
    label 'compute_resources__sratoolkit'

    input:
        tuple val(sraId), val(sampleId)
    
    output:
        tuple val(sraId), file("${sraId}*.fastq")
    
    script:
        if(sampleId == null || sampleId.length() < 1) {
            throw new Exception("DOWNLOAD_FASTQS_FROM_SRA_ACC_ID: Sample ID is empty.")
        }
        """
        SRA_FILE_LOCK=./ncbi/public/sra/${sraId}.sra.lock
        if [[ -f "\${SRA_FILE_LOCK}" ]]; then
            echo "SRA file lock found for ${sraId}. Removing file lock..."
            rm \${SRA_FILE_LOCK}
        fi
        # Fetch SRA file
        prefetch \
           -v \
           -p 1 \
           ${params.tools.sratoolkit?.maxSize ? '--max-size '+ params.tools.sratoolkit.maxSize: ''} \
           ${sraId}
        # Convert SRA file to FASTQ files
        fasterq-dump \
           -S \
           -v \
           -p \
           -e ${task.cpus} \
           ${params.tools.sratoolkit?.includeTechnicalReads ? '--include-technical' : ''} \
           -O . \
           ${sraId}
        """

}
