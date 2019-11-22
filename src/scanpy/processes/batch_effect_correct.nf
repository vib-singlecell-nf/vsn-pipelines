nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  binDir = ""
}

process SC__SCANPY__BATCH_EFFECT_CORRECTION {

  container params.sc.scanpy.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=6gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true
  
  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.SC__SCANPY__BATCH_EFFECT_CORRECTION.${params.off}")
  script:
    """
    ${binDir}aggregate/sc_batch_effect_correction.py \
        ${(params.containsKey('batchEffectCorrectionMethod')) ? '--method ' + params.batchEffectCorrectionMethod : ''} \
        --output-file "${id}.SC__SCANPY__BATCH_EFFECT_CORRECTION.${params.off}" \
        ${(params.containsKey('key')) ? '--key ' + params.key : ''} \
        ${(params.containsKey('batchKey')) ? '--batch-key ' + params.batchKey : ''} \
        ${(params.containsKey('nPcs')) ? '--n-pcs ' + params.nPcs : ''} \
        ${(params.containsKey('k')) ? '--k' + params.k : ''} \
        ${(params.containsKey('varIndex')) ? '--var-index ' + params.varIndex : ''} \
        ${(params.containsKey('varSubset')) ? '--var-subset ' + params.varSubset : ''} \
        ${(params.containsKey('nJobs')) ? '--n-jobs ' + params.nJobs : ''} \
        ${(params.containsKey('neighborsWithinBatch')) ? '--neighbors-within-batch ' + params.neighborsWithinBatch : ''} \
        ${(params.containsKey('trim')) ? '--trim ' + params.trim : ''} \
        $f
    """
}
