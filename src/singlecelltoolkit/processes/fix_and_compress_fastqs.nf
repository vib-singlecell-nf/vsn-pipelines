nextflow.enable.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/sratoolkit/bin/"
} else {
    binDir = ""
}

toolParams = params.tools.sratoolkit

process FIX_AND_COMPRESS_SRA_FASTQS {

    container toolParams.container
    publishDir "${params.global.outdir}/data/raw/fastqs_fixed_and_compressed", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sraId), file("${sraId}_*.fastq")
    
    output:
        tuple val(sraId), file("${sraId}_*.fastq.gz")
    
    script:
        """
        # Fixing the FASTQ files is required for future pre-processing (e.g.: scATAC-seq pipelines) because fasterq-dump does not have the -F option as fastq-dump do to keep original sequence names.
        # Fix the FASTQ files and compress them
        export compress_fastq_threads="${task.cpus}"
        NUM_FASTQ_FILES=\$(ls ./*.fastq | wc -l)
        echo "Fixing and compressing \${NUM_FASTQ_FILES} FASTQ files in parallel with \${compress_fastq_threads} compression threads for each task..."
        echo *.fastq | tr ' ' '\n' | xargs -P "\${NUM_FASTQ_FILES}" -n 1 -I {} fix_sra_fastq.sh "{}" "{}.gz" pigz
        echo "Removing all uncompressed FASTQ files"
        for FASTQ in *.fastq; do
           echo "Removing uncompressed FASTQ file \${FASTQ}..."
           rm "\$(readlink -f \${FASTQ})"
        done
        echo "Done."
        """

}
