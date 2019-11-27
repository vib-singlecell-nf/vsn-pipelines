nextflow.preview.dsl=2

/**
 * Preprocess + FastQC
 */
process FASTP__CLEAN_AND_FASTQC {

    container params.fastp.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    set val(sample), path(reads)
    
    output:
    tuple file('*_R{1,2}.clean.fastq.gz'), emit: fastq
    tuple file('*_fastp.{json,html}'), emit: report
    
    script:
    processParams = params.fastp
    """
    fastp --thread ${processParams.thread} \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${sample}_R1.clean.fastq.gz \
        -O ${sample}_R2.clean.fastq.gz \
        --length_required ${processParams.clean_and_fastqc.length_required} \
        --adapter_fasta ${processParams.clean_and_fastqc.adapter_fasta} \
        -j ${sample}_fastp.json \
        -h ${sample}_fastp.html
    """
}
