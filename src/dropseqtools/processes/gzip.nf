nextflow.preview.dsl=2

process GZIP {

    container params.sc.dropseqtools.container
    publishDir "${params.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(f)
    output:
        tuple val(sample), file("*.unaligned_tagged_polyA_filtered.fastq.gz"), emit: fastq_gz
    script:
        """
        gzip --force ${f}
        """
}