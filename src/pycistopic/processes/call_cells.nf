nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/pycistopic/bin/" : ""

toolParams = params.tools.pycistopic

process SC__PYCISTOPIC__CALL_CELLS {

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
              path(output_pdf)

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.call_cells)
        processParams = sampleParams.local
        selected_barcodes = "${sampleId}__selected_barcodes.txt"
        output_metadata_pkl = "${sampleId}__metadata_with_calls.pickle"
        output_pdf = "${sampleId}__fragments_qc.pdf"
        """
        ${binDir}call_cells.py \
            --sampleId ${sampleId} \
            --metadata_pkl ${metadata_pkl} \
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

