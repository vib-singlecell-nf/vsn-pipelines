nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

process SC__SCANPY__BATCH_EFFECT_CORRECTION {

  	container params.getToolParams("scanpy").container
  	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

  	input:
    	tuple \
			val(sampleId), \
			path(f), \
			val(stashedParams)

  	output:
    	tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCANPY__BATCH_EFFECT_CORRECTION.${processParams.off}"), \
			val(stashedParams)

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("scanpy").batch_effect_correct)
		processParams = sampleParams.local
		"""
		${binDir}aggregate/sc_batch_effect_correction.py \
			${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
			--output-file "${sampleId}.SC__SCANPY__BATCH_EFFECT_CORRECTION.${processParams.off}" \
			${(processParams.containsKey('key')) ? '--key ' + processParams.key : ''} \
			${(processParams.containsKey('batchKey')) ? '--batch-key ' + processParams.batchKey : ''} \
			${(processParams.containsKey('nPcs')) ? '--n-pcs ' + processParams.nPcs : ''} \
			${(processParams.containsKey('k')) ? '--k ' + processParams.k : ''} \
			${(processParams.containsKey('varIndex')) ? '--var-index ' + processParams.varIndex : ''} \
			${(processParams.containsKey('varSubset')) ? '--var-subset ' + processParams.varSubset : ''} \
            --n-jobs ${task.cpus} \
			${(processParams.containsKey('neighborsWithinBatch')) ? '--neighbors-within-batch ' + processParams.neighborsWithinBatch : ''} \
			${(processParams.containsKey('trim')) ? '--trim ' + processParams.trim : ''} \
			$f
		"""

}
