nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/popscle/bin/" : ""

process SC__POPSCLE__DSC_PILEUP {

    container params.sc.popscle.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'

    input:
        tuple val(sampleId), path(f)

    // output:
    //     tuple val(sampleId), path("${sampleId}.SC__POPSCLE__DSC_PILEUP.h5ad")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.popscle.dsc_pileup)
		processParams = sampleParams.local

        """
        popscle dsc-pileup \
            --sam ${f}/possorted_genome_bam.bam \
            --vcf ${processParams.vcf} \
            --out ${sampleId}_dsc-pileup
        """
}

