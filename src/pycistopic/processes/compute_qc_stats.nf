nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/pycistopic/bin/" : ""

toolParams = params.tools.pycistopic
processParams = params.tools.pycistopic.compute_qc_stats

process rename_fragments {

    container toolParams.container
    label 'compute_resources__minimal'

    input:
        tuple val(sampleId),
              path(f)
    output:
        tuple val(sampleId),
              path("${sampleId}_${f}*")

    script:
        """
        ln -s ${f[0]} ${sampleId}_${f[0]}
        ln -s ${f[1]} ${sampleId}_${f[1]}
        """

}


process PYCISTOPIC__COMPUTE_QC_STATS {

    publishDir "${params.global.outdir}/data/pycistopic/qc/", mode: params.utils.publish.mode
    container toolParams.container
    label 'compute_resources__pycisTopic'

    input:
        val(input)
        path(biomart_annot)
        path(fragments)
        path(peaks)

    output:
        tuple path('metadata/*.metadata.pkl'),
              path('profile_data/*.profile_data.pkl')

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
            --output_metadata_dir metadata \
            --output_profile_data_dir profile_data
        """
}

