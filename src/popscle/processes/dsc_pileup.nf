nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/popscle/bin/" : ""

toolParams = params.tools.popscle

process SC__POPSCLE__DSC_PILEUP {

    container params.tools.popscle.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sampleId), path(f)
        file vcf

    output:
        tuple val(sampleId), path("${sampleId}_dsc-pileup*.gz")

    script:
        """
        popscle dsc-pileup \
            --sam ${f} \
            ${toolParams?.barcode_tag ? '--tag-group ' +  toolParams.barcode_tag : ''} \
            --vcf ${vcf} \
            --out ${sampleId}_dsc-pileup
        """
}

process SC__POPSCLE__PREFILTER_DSC_PILEUP {

    container params.tools.popscle.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId),
              path(bam),
              path(barcodes)
        file vcf

    output:
        tuple val(sampleId), path("${sampleId}_filtered_possorted_genome_bam.bam")

    script:
        """
        filter_bam_file_for_popscle_dsc_pileup.sh \
            ${bam} \
            ${barcodes} \
            ${vcf} \
            ${sampleId}_filtered_possorted_genome_bam.bam \
            ${toolParams?.barcode_tag ? toolParams.barcode_tag : ''}
        """
}

