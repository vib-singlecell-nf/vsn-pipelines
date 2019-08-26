nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

process SC__SCANPY__BATCH_EFFECT_CORRECTION {
  input:
    file(f)
  output:
    file "${params.project_name}.SC__SCANPY__BATCH_EFFECT_CORRECTION.${params.off}" //#"${getBaseName(f).get()}.SC__SCANPY__BATCH_EFFECT_CORRECTION.${params.off}" 
  script:
    """
    python ../../../src/scripts/scanpy/aggregate/sc_batch_effect_correction.py \
        ${(params.containsKey('batchEffectCorrectionMethod')) ? '--method ' + params.batchEffectCorrectionMethod : ''} \
        --output-file "${params.project_name}.SC__SCANPY__BATCH_EFFECT_CORRECTION.${params.off}" \
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