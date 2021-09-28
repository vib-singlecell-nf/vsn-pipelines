nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/pycistopic/bin/" : ""

toolParams = params.tools.pycistopic
processParams = params.tools.pycistopic.compute_qc_stats

process PYCISTOPIC__PLOT_QC_STATS {

    container toolParams.container
    label 'compute_resources__default'

    input:
        tuple val(sampleId),
              path(output_metadata),
              path(output_metadata_pkl),
              path(output_profile_data_pkl)

    output:
        tuple val(sampleId),
              path(output_pdf)

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        output_metadata = "${sampleId}_metadata.tsv.gz"
        output_pdf = "${sampleId}_qc_sample_metrics.pdf"
        output_metadata_pkl = "${sampleId}_metadata.pickle"
        output_profile_data_pkl = "${sampleId}_profile_data.pickle"
        """
        ${binDir}plot_qc_stats.py \
            --sampleId ${sampleId} \
            --profile_data_pkl ${output_profile_data_pkl} \
            --output_pdf ${output_pdf}
        """
}

