nextflow.enable.dsl=2

/**
 * Preprocess + FastQC
 */
process FASTP__CLEAN_AND_FASTQC {

    container params.tools.fastp.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sampleId), path(reads)
    
    output:
        tuple val(sampleId), path('*_R{1,2}.clean.fastq.gz'), emit: fastq
        tuple val(sampleId), path('*_fastp.{json,html}'), emit: report
    
    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.fastp)
		processParams = sampleParams.local

        """
        case "${reads[0]}" in
            *.fastq.gz) EXTENSION='fastq.gz' ;;
            *.fastq)     EXTENSION='fastq' ;;
            *) ;;
        esac

        if [ `ls -1 ${sampleId}_L00[0-9]_R1_*.\$EXTENSION 2>/dev/null | wc -l ` -gt 1 ]; then
            cat ${sampleId}_L00[0-9]_R1_*.\$EXTENSION > ${sampleId}_R1.\$EXTENSION
        else
            mv ${reads}[0] ${sampleId}_R1.\$EXTENSION
        fi

        if [ `ls -1 ${sampleId}_L00[0-9]_R2_*.\$EXTENSION 2>/dev/null | wc -l ` -gt 1 ]; then
            cat ${sampleId}_L00[0-9]_R2_*.\$EXTENSION > ${sampleId}_R2.\$EXTENSION
        else
            mv ${reads}[0] ${sampleId}_R2.\$EXTENSION
        fi

        fastp \
            -i ${sampleId}_R1.\$EXTENSION \
            -I ${sampleId}_R2.\$EXTENSION \
            -o ${sampleId}_R1.clean.fastq.gz \
            -O ${sampleId}_R2.clean.fastq.gz \
            --thread ${task.cpus} \
            --length_required ${processParams.clean_and_fastqc.length_required} \
            --adapter_fasta ${processParams.clean_and_fastqc.adapter_fasta} \
            -j ${sampleId}_fastp.json \
            -h ${sampleId}_fastp.html
        """
}
