nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

params.adjustmentMethod = 'linear_regression'
params.normalizationVariablesToRegressOut = ['n_counts','percent_mito']
normalizationVariablesToRegressOutAsArguments = params.normalizationVariablesToRegressOut.collect({ '--variable-to-regress-out' + ' ' + it }).join(' ')

process SC__SCANPY__ADJUSTMENT {
  input:
    file(f)
  output:
    file "${getBaseName(f).get()}.SC__SCANPY__ADJUSTMENT.${params.off}" 
  script:
    """
    python ../../../src/scripts/scanpy/adjust/sc_adjustment.py \
         --method $params.adjustmentMethod \
         $normalizationVariablesToRegressOutAsArguments \
         $f \
         "${getBaseName(f).get()}.SC__SCANPY__ADJUSTMENT.${params.off}" 
    """
}
