nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

process SC__SCANPY__REGRESS_OUT {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

	input:
		tuple val(sampleId), path(f)

	output:
		tuple val(sampleId), path("${sampleId}.SC__SCANPY__REGRESS_OUT.${params.off}")

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.regress_out)
		processParams = sampleParams.local
		variablesToRegressOutAsArguments = processParams.variablesToRegressOut.collect({ '--variable-to-regress-out' + ' ' + it }).join(' ')
		"""
		${binDir}adjust/sc_regress_out.py \
			${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
			${(processParams.containsKey('variablesToRegressOut')) ? variablesToRegressOutAsArguments : ''} \
			$f \
			"${sampleId}.SC__SCANPY__REGRESS_OUT.${params.off}" 
		"""

}
