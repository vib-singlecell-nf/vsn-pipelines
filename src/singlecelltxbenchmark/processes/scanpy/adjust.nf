nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

params.normalizationVariablesToRegressOut = ['n_counts','percent_mito']
normalizationVariablesToRegressOutAsArguments = params.normalizationVariablesToRegressOut.collect({ '--variable-to-regress-out' + ' ' + it }).join(' ')

process SC__SCANPY__ADJUSTMENT {
  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__ADJUSTMENT.${params.off}" 
  script:
    """
    python ../../../src/singlecelltxbenchmark/scripts/scanpy/adjust/sc_adjustment.py \
         ${(params.containsKey('adjustmentMethod')) ? '--method ' + params.adjustmentMethod : ''} \
         ${(params.containsKey('normalizationVariablesToRegressOut')) ? normalizationVariablesToRegressOutAsArguments : ''} \
         $f \
         "${getBaseName(f)}.SC__SCANPY__ADJUSTMENT.${params.off}" 
    """
}
