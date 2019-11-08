nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process SC__SCENIC__AGGR_MULTI_RUNS_FEATURES {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/multi_runs_cistarget/", mode: 'link', overwrite: true
    // In the case the chunking method is not used, this process requires a large amount of memory especially for big datasets
    // This process is quite slow (could take more than 1h for big datasets, so keep 24h for now)
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=6gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    file f
    val type

    output:
    file "multi_runs_features_${type}.${output_format_ext}${compression_ext}"

    script:
    output_format = params.sc.scenic.aggregate_features.output_format
    output_format_ext = output_format
    if(output_format == 'pickle') {
      output_format_ext = 'pkl'
    }
    compression = params.sc.scenic.aggregate_features.compression
    compression_ext = ''
    if(compression == 'gzip') {
      compression_ext = '.gz'
    }
    """
    ${binDir}aggregate_SCENIC_multi_runs_features.py \
        ${f} \
        --output "multi_runs_features_${type}.${output_format_ext}${compression_ext}" \
        --output-format ${output_format} \
        --use-chunking ${params.sc.scenic.aggregate_features.use_chunking}
    """
}
