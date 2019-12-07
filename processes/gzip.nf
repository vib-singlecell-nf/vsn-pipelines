nextflow.preview.dsl=2

process GZIP {

    container params.sc.dropseqtools.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sample), path(f)

    output:
        tuple val(sample), path("*.unaligned_tagged_polyA_filtered.fastq.gz"), emit: fastq_gz

    script:
        """
        gzip --force ${f}
        """

}