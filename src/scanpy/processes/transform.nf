nextflow.preview.dsl=2

if(!params.containsKey("test")) {
	binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
	binDir = ""
}

process SC__SCANPY__NORMALIZATION {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

	input:
		tuple val(sampleId), path(f)

	output:
		tuple val(sampleId), path("${sampleId}.SC__SCANPY__NORMALIZATION.${processParams.off}")

	script:
		processParams = params.sc.scanpy.normalization
		"""
		${binDir}transform/sc_normalization.py \
			${(processParams.containsKey('normalizationMethod')) ? '--method ' + processParams.normalizationMethod : ''} \
			${(processParams.containsKey('countsPerCellAfter')) ? '--counts-per-cell-after ' + processParams.countsPerCellAfter : ''} \
			$f \
			"${sampleId}.SC__SCANPY__NORMALIZATION.${processParams.off}"
		"""

}

process SC__SCANPY__DATA_TRANSFORMATION {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

	input:
		tuple val(sampleId), path(f)
	
	output:
		tuple val(sampleId), path("${sampleId}.SC__SCANPY__DATA_TRANSFORMATION.${processParams.off}")
	
	script:
		processParams = params.sc.scanpy.data_transformation
		"""
		${binDir}transform/sc_data_transformation.py \
			${(processParams.containsKey('dataTransformationMethod')) ? '--method ' + processParams.dataTransformationMethod : ''} \
			$f \
			"${sampleId}.SC__SCANPY__DATA_TRANSFORMATION.${processParams.off}"
		"""

}

process SC__SCANPY__FEATURE_SCALING {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

	input:
		tuple val(sampleId), path(f)
	
	output:
		tuple val(sampleId), path("${sampleId}.SC__SCANPY__FEATURE_SCALING.${processParams.off}")
	
	script:
		processParams = params.sc.scanpy.feature_scaling
		"""
		${binDir}transform/sc_feature_scaling.py \
			${(processParams.containsKey('featureScalingMthod')) ? '--method ' + processParams.featureScalingMthod : ''} \
			${(processParams.containsKey('featureScalingMaxSD')) ? '--max-sd ' + processParams.featureScalingMaxSD : ''} \
			$f \
			"${sampleId}.SC__SCANPY__FEATURE_SCALING.${processParams.off}"
		"""

}
