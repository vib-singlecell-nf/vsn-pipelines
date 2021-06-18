nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/pycistopic/bin/" : ""

toolParams = params.tools.pycistopic
processParams = params.tools.pycistopic.compute_qc_stats

process PYCISTOPIC__COMPUTE_QC_STATS {

    publishDir "${params.global.outdir}/data/pycistopic/qc/", mode: params.utils.publish.mode
    container toolParams.container
    label 'compute_resources__pycisTopic'

    input:
        val(input)
        path(biomart_annot)

    output:
        tuple path("project_metadata.pickle"),
              path("project_profile_data.pickle")

    script:
        """
        export NUMEXPR_MAX_THREADS=1
        export OMP_NUM_THREADS=1
        ${binDir}compute_qc_stats.py \
            ${"--input_files "+input.join(" --input_files ")} \
            --n_frag ${processParams.n_frag} \
            --tss_flank_window ${processParams.tss_flank_window} \
            --tss_window ${processParams.tss_window} \
            --tss_minimum_signal_window ${processParams.tss_minimum_signal_window} \
            --tss_rolling_window ${processParams.tss_rolling_window} \
            --min_norm ${processParams.min_norm} \
            --threads ${task.cpus} \
            --biomart_annot_pkl ${biomart_annot} \
            --output_metadata_pkl project_metadata.pickle \
            --output_profile_data_pkl project_profile_data.pickle
        """
}

