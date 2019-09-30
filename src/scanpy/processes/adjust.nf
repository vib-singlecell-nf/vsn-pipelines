nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

params.normalizationVariablesToRegressOut = ['n_counts','percent_mito']
normalizationVariablesToRegressOutAsArguments = params.normalizationVariablesToRegressOut.collect({ '--variable-to-regress-out' + ' ' + it }).join(' ')

process SC__SCANPY__ADJUSTMENT {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/data", mode: 'symlink'

  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__ADJUSTMENT.${params.off}" 
  script:
    """
    ${workflow.projectDir}/src/scanpy/bin/adjust/sc_adjustment.py \
         ${(params.containsKey('adjustmentMethod')) ? '--method ' + params.adjustmentMethod : ''} \
         ${(params.containsKey('normalizationVariablesToRegressOut')) ? normalizationVariablesToRegressOutAsArguments : ''} \
         $f \
         "${getBaseName(f)}.SC__SCANPY__ADJUSTMENT.${params.off}" 
    """
}
