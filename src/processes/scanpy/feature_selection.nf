nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

process SC__SCANPY__FEATURE_SELECTION {
  input:
    file(f)
  output:
    file "${getBaseName(f).get()}.SC__SCANPY__FEATURE_SELECTION.${params.off}"
  script:
    """
    python ../../../src/scripts/scanpy/feature_selection/sc_select_variable_genes.py \
        --method $params.featureSelectionMethod \
        ${(params.containsKey('featureSelectionMinMean')) ? '--min-mean ' + params.featureSelectionMinMean : ''} \
        ${(params.containsKey('featureSelectionMaxMean')) ? '--max-mean ' + params.featureSelectionMaxMean : ''} \
        ${(params.containsKey('featureSelectionMinDisp')) ? '--min-disp ' + params.featureSelectionMinDisp : ''} \
        ${(params.containsKey('featureSelectionMaxDisp')) ? '--max-disp ' + params.featureSelectionMaxDisp : ''} \
        $f \
        "${getBaseName(f).get()}.SC__SCANPY__FEATURE_SELECTION.${params.off}"
    """
}