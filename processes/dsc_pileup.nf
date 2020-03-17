nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/popscle/bin/" : ""

process SC__POPSCLE__DSC_PILEUP {

    container params.sc.popscle.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'

    input:
        tuple val(sampleId), path(f)
        file vcf

    output:
        tuple val(sampleId), path("${sampleId}_dsc-pileup*.gz")

    script:
        """
        popscle dsc-pileup \
            --sam ${f} \
            --vcf ${vcf} \
            --out ${sampleId}_dsc-pileup
        """
}

process SC__POPSCLE__PREFILTER_DSC_PILEUP {

    container params.sc.popscle.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'

    input:
        tuple val(sampleId), path(f)
        file vcf

    output:
        tuple val(sampleId), path("${sampleId}_filtered_possorted_genome_bam.bam")

    script:
        """
        ${binDir}/filter_bam_file_for_popscle_dsc_pileup.sh \
            ${f}/possorted_genome_bam.bam \
            ${f}/filtered_*_bc_matrix/barcodes.tsv* \
            ${vcf} \
            ${sampleId}_filtered_possorted_genome_bam.bam
        """
}
