nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/pycistopic/bin/" : ""

toolParams = params.tools.pycistopic
//processParams = params.tools.pycistopic.barcode_level_statistics

process PYCISTOPIC__BARCODE_LEVEL_STATISTICS {

    publishDir "${params.global.outdir}/intermediate/pycistopic/qc/", mode: 'symlink'
    container toolParams.container
    label 'compute_resources__default','compute_resources__24hqueue'

    input:
        tuple val(sampleId),
              path(metadata),
              path(metadata_pkl),
              path(profile_data_pkl)

    output:
        tuple val(sampleId),
              path(selected_barcodes),
              path(output_pdf_ff),
              path(output_pdf_tf),
              path(output_pdf_df)

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.barcode_level_statistics)
        processParams = sampleParams.local
        selected_barcodes = "${sampleId}__selected_barcodes.txt"
        output_pdf_ff = "${sampleId}__FRIP-vs-nFrag.pdf"
        output_pdf_tf = "${sampleId}__TSS-vs-nFrag.pdf"
        output_pdf_df = "${sampleId}__duprate-vs-nFrag.pdf"
        """
        export NUMEXPR_MAX_THREADS=${task.cpus}
        ${binDir}barcode_level_statistics.py \
            --sampleId ${sampleId} \
            --metadata_pkl ${metadata_pkl} \
            --selected_barcodes ${selected_barcodes} \
            ${processParams?.filter_frags_lower    ? '--filter_frags_lower '    + processParams?.filter_frags_lower    : ''} \
            ${processParams?.filter_frags_upper    ? '--filter_frags_upper '    + processParams?.filter_frags_upper    : ''} \
            ${processParams?.filter_tss_lower      ? '--filter_tss_lower '      + processParams?.filter_tss_lower      : ''} \
            ${processParams?.filter_tss_upper      ? '--filter_tss_upper '      + processParams?.filter_tss_upper      : ''} \
            ${processParams?.filter_frip_lower     ? '--filter_frip_lower '     + processParams?.filter_frip_lower     : ''} \
            ${processParams?.filter_frip_upper     ? '--filter_frip_upper '     + processParams?.filter_frip_upper     : ''} \
            ${processParams?.filter_dup_rate_lower ? '--filter_dup_rate_lower ' + processParams?.filter_dup_rate_lower : ''} \
            ${processParams?.filter_dup_rate_upper ? '--filter_dup_rate_upper ' + processParams?.filter_dup_rate_upper : ''}
        """
}

