nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  binDir = ""
}

process SC__SCANPY__ADJUSTMENT {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

	input:
	tuple val(id), file(f)

	output:
	tuple val(id), file("${id}.SC__SCANPY__ADJUSTMENT.${params.off}")

	script:
	processParams = params.sc.scanpy.data_adjustment
	normalizationVariablesToRegressOutAsArguments = processParams.normalizationVariablesToRegressOut.collect({ '--variable-to-regress-out' + ' ' + it }).join(' ')
	"""
	${binDir}adjust/sc_adjustment.py \
		${(processParams.containsKey('adjustmentMethod')) ? '--method ' + processParams.adjustmentMethod : ''} \
		${(processParams.containsKey('normalizationVariablesToRegressOut')) ? normalizationVariablesToRegressOutAsArguments : ''} \
		$f \
		"${id}.SC__SCANPY__ADJUSTMENT.${params.off}" 
	"""

}
