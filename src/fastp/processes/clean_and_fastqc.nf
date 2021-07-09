nextflow.enable.dsl=2

/**
 * Preprocess + FastQC
 */
process FASTP__CLEAN_AND_FASTQC {

    container params.tools.fastp.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        set val(sample), path(reads)
    
    output:
        tuple file('*_R{1,2}.clean.fastq.gz'), emit: fastq
        tuple file('*_fastp.{json,html}'), emit: report
    
    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.fastp)
		processParams = sampleParams.local
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
