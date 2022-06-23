nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.bap

process BAP__BARCODE_MULTIPLET_PIPELINE {

    container toolParams.container
    publishDir "${params.global.outdir}/data/bap", mode: params.utils.publish.mode
    label 'compute_resources__bap_barcode_multiplet_pipeline'

    input:
        tuple val(sampleId),
              path("input.bam"),
              path("input.bam.bai")

    output:
        tuple val(sampleId),
              path("${sampleId}/final/${sampleId}.bap.bam"),
              path("${sampleId}/final/${sampleId}.bap.bam.bai"),
              path("${sampleId}/final/*"),
              path("${sampleId}/knee/*"),
              path("${sampleId}/logs/*")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.barcode_multiplet)
        processParams = sampleParams.local
        """
        export OMP_THREAD_LIMIT=8
        bap2 bam \
            --input input.bam \
            --output ${sampleId} \
            --name ${sampleId} \
            --ncores ${task.cpus} \
            --drop-tag ${processParams.drop_tag} \
            --bead-tag ${processParams.bead_tag} \
            --minimum-barcode-fragments ${processParams.minimum_barcode_fragments} \
            ${processParams?.barcode_whitelist ? '--barcode-whitelist ' + processParams.barcode_whitelist : ''} \
            --minimum-jaccard-index ${processParams.minimum_jaccard_index} \
            --nc-threshold ${processParams.nc_threshold} \
            --mapq ${processParams.mapq} \
            --max-insert ${processParams.max_insert} \
            ${processParams?.reference_genome ? '--reference-genome ' + processParams.reference_genome : ''} \
            ${processParams?.bedtools_genome ? '--bedtools-genome ' + processParams.bedtools_genome : ''} \
            ${processParams?.blacklist_file ? '--blacklist-file ' + processParams.blacklist_file : ''} \
            ${processParams?.tss_file ? '--tss-file ' + processParams.tss_file : ''} \
            --mito-chromosome ${processParams.mito_chromosome}
        """
}

