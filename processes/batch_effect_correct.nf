nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  	binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  	binDir = ""
}

process SC__SCANPY__BATCH_EFFECT_CORRECTION {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=6gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
  
  	input:
    tuple val(sampleId), path(f)

  	output:
    tuple val(sampleId), path("${sampleId}.SC__SCANPY__BATCH_EFFECT_CORRECTION.${processParams.off}")
  	
	script:
    processParams = params.sc.scanpy.batch_effect_correct
    """
    ${binDir}aggregate/sc_batch_effect_correction.py \
        ${(processParams.containsKey('batchEffectCorrectionMethod')) ? '--method ' + processParams.batchEffectCorrectionMethod : ''} \
        --output-file "${sampleId}.SC__SCANPY__BATCH_EFFECT_CORRECTION.${processParams.off}" \
        ${(processParams.containsKey('key')) ? '--key ' + processParams.key : ''} \
        ${(processParams.containsKey('batchKey')) ? '--batch-key ' + processParams.batchKey : ''} \
        ${(processParams.containsKey('nPcs')) ? '--n-pcs ' + processParams.nPcs : ''} \
        ${(processParams.containsKey('k')) ? '--k' + processParams.k : ''} \
        ${(processParams.containsKey('varIndex')) ? '--var-index ' + processParams.varIndex : ''} \
        ${(processParams.containsKey('varSubset')) ? '--var-subset ' + processParams.varSubset : ''} \
        ${(processParams.containsKey('nJobs')) ? '--n-jobs ' + processParams.nJobs : ''} \
        ${(processParams.containsKey('neighborsWithinBatch')) ? '--neighbors-within-batch ' + processParams.neighborsWithinBatch : ''} \
        ${(processParams.containsKey('trim')) ? '--trim ' + processParams.trim : ''} \
        $f
    """
}
