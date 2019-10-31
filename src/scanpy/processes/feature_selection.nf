nextflow.preview.dsl=2

include getBaseName from '../../utils/processes/files.nf'

process SC__SCANPY__FEATURE_SELECTION {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__FEATURE_SELECTION.${params.off}"
  script:
    """
    ${workflow.projectDir}/src/scanpy/bin/feature_selection/sc_select_variable_genes.py \
        --method $params.featureSelectionMethod \
        ${(params.containsKey('featureSelectionMinMean')) ? '--min-mean ' + params.featureSelectionMinMean : ''} \
        ${(params.containsKey('featureSelectionMaxMean')) ? '--max-mean ' + params.featureSelectionMaxMean : ''} \
        ${(params.containsKey('featureSelectionMinDisp')) ? '--min-disp ' + params.featureSelectionMinDisp : ''} \
        ${(params.containsKey('featureSelectionMaxDisp')) ? '--max-disp ' + params.featureSelectionMaxDisp : ''} \
        $f \
        "${getBaseName(f)}.SC__SCANPY__FEATURE_SELECTION.${params.off}"
    """
}

