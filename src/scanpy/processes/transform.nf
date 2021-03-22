nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SCANPY__NORMALIZATION {

	container params.sc.scanpy.container
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

	input:
		tuple val(sampleId), path(f)

	output:
		tuple val(sampleId), path("${sampleId}.SC__SCANPY__NORMALIZATION.${processParams.off}")

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.normalization)
		processParams = sampleParams.local
		"""
		${binDir}/transform/sc_normalization.py \
			${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
			${(processParams.containsKey('targetSum')) ? '--target-sum ' + processParams.targetSum : ''} \
			$f \
			"${sampleId}.SC__SCANPY__NORMALIZATION.${processParams.off}"
		"""

}

process SC__SCANPY__DATA_TRANSFORMATION {

	container params.sc.scanpy.container
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

	input:
		tuple val(sampleId), path(f)
	
	output:
		tuple val(sampleId), path("${sampleId}.SC__SCANPY__DATA_TRANSFORMATION.${processParams.off}")
	
	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.data_transformation)
		processParams = sampleParams.local
		"""
		${binDir}/transform/sc_data_transformation.py \
			${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
			$f \
			"${sampleId}.SC__SCANPY__DATA_TRANSFORMATION.${processParams.off}"
		"""

}

process SC__SCANPY__FEATURE_SCALING {

	container params.sc.scanpy.container
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

	input:
		tuple \
			val(sampleId), \
			path(f)
	
	output:
		tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCANPY__FEATURE_SCALING.${processParams.off}")
	
	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.feature_scaling)
		processParams = sampleParams.local
		"""
		${binDir}/transform/sc_feature_scaling.py \
			${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
			${(processParams.containsKey('maxSD')) ? '--max-sd ' + processParams.maxSD : ''} \
			$f \
			"${sampleId}.SC__SCANPY__FEATURE_SCALING.${processParams.off}"
		"""

}
