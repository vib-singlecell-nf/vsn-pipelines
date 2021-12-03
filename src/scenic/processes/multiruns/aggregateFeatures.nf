nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

def toolParams = params.tools.scenic
def processParams = toolParams.aggregate_features

process AGGR_MULTI_RUNS_FEATURES {

    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/multi_runs_cistarget/", mode: 'link', overwrite: true
    // In the case the chunking method is not used, this process requires a large amount of memory especially for big datasets
    // This process is quite slow (could take more than 1h for big datasets, so keep 24h for now)
    label 'compute_resources__scenic_multiruns'

    input:
		tuple val(sampleId), path(f)
		val type

    output:
    	tuple val(sampleId), path("multi_runs_features_${type}.${output_format_ext}${compression_ext}")

    script:
		output_format = processParams.output_format
		output_format_ext = output_format
		if(output_format == 'pickle') {
		output_format_ext = 'pkl'
		}
		compression = processParams.compression
		compression_ext = ''
		if(compression == 'gzip') {
		compression_ext = '.gz'
		}
		"""
		${binDir}aggregate_multi_runs_features.py \
			${f} \
			--output "multi_runs_features_${type}.${output_format_ext}${compression_ext}" \
			--output-format ${output_format} \
			--use-chunking ${processParams.use_chunking}
		"""

}
