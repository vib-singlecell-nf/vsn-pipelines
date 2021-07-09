nextflow.enable.dsl=2

process GZIP {

    container params.tools.dropseqtools.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sample), path(f)

    output:
        tuple val(sample), path("*.unaligned_tagged_polyA_filtered.fastq.gz"), emit: fastq_gz

    script:
        """
        gzip --force ${f}
        """

}
