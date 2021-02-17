nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

process SC__SCANPY__REGRESS_OUT {

	container params.tools.scanpy.container
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

	input:
		tuple \
			val(sampleId), \
			path(f)

	output:
		tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCANPY__REGRESS_OUT.${processParams.off}")

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.tools.scanpy.regress_out)
		processParams = sampleParams.local
		variablesToRegressOutAsArguments = processParams.variablesToRegressOut.collect({ '--variable-to-regress-out' + ' ' + it }).join(' ')
		"""
		export MKL_NUM_THREADS=1
		export NUMEXPR_NUM_THREADS=1
		export OMP_NUM_THREADS=1
		${binDir}adjust/sc_regress_out.py \
			${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
			${(processParams.containsKey('variablesToRegressOut')) ? variablesToRegressOutAsArguments : ''} \
			--n-jobs ${task.cpus} \
			$f \
			"${sampleId}.SC__SCANPY__REGRESS_OUT.${processParams.off}"
		"""

}
