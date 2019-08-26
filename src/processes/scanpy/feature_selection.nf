nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

params.featureSelectionMethod = 'mean_disp_plot'
params.featureSelectionMinMean = -1
params.featureSelectionMaxMean = -1
params.featureSelectionMinDisp = -1
params.featureSelectionMaxDisp = -1

process SC__SCANPY__FEATURE_SELECTION {
  input:
    file(f)
  output:
    file "${getBaseName(f).get()}.SC__SCANPY__FEATURE_SELECTION.${params.off}"
  script:
    """
    python ../../../src/scripts/scanpy/feature_selection/sc_select_variable_genes.py \
        --method $params.featureSelectionMethod \
        --min-mean $params.featureSelectionMinMean \
        --max-mean $params.featureSelectionMaxMean \
        --min-disp $params.featureSelectionMinDisp \
        --max-disp $params.featureSelectionMaxDisp \
       $f \
      "${getBaseName(f).get()}.SC__SCANPY__FEATURE_SELECTION.${params.off}"
    """
}