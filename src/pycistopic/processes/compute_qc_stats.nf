nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/pycistopic/bin/" : ""

toolParams = params.tools.pycistopic
processParams = params.tools.pycistopic.compute_qc_stats

process SC__PYCISTOPIC__COMPUTE_QC_STATS {

    container toolParams.container
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId),
              path(fragments),
              path(fragments_index),
              val(peaks)

    output:
        tuple val(sampleId),
              path(output_metadata),
              path(output_metadata_pkl),
              path(output_profile_data_pkl)

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        output_metadata = "${sampleId}_metadata.tsv.gz"
        output_metadata_pkl = "${sampleId}_metadata.pickle"
        output_profile_data_pkl = "${sampleId}_profile_data.pickle"
        """
        export NUMEXPR_MAX_THREADS=${task.cpus}
        ${binDir}compute_qc_stats.py \
            --sampleId ${sampleId} \
            --fragments ${fragments} \
            --regions ${peaks} \
            --n_frag ${processParams.n_frag} \
            --threads ${task.cpus} \
            --output_metadata ${output_metadata} \
            --output_metadata_pkl ${output_metadata_pkl} \
            --output_profile_data_pkl ${output_profile_data_pkl}
        """
}

