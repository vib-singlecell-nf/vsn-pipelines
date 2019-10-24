nextflow.preview.dsl=2

process PICARD__BAM_TO_FASTQ {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_polyA_filtered.fastq'), emit: fastq
    script:
        """
        module load Java/1.8.0_192
        java -Djava.io.tmpdir=$DWMAX/tmp -jar \
            /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/picard/2.20.6/bin/picard.jar \
                SamToFastq \
			        INPUT=${bam} \
			        FASTQ=${sample}.unaligned_tagged_polyA_filtered.fastq
        """    
}