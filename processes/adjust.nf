nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  binDir = ""
}

params.normalizationVariablesToRegressOut = ['n_counts','percent_mito']
normalizationVariablesToRegressOutAsArguments = params.normalizationVariablesToRegressOut.collect({ '--variable-to-regress-out' + ' ' + it }).join(' ')

process SC__SCANPY__ADJUSTMENT {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.SC__SCANPY__ADJUSTMENT.${params.off}")
  script:
    """
    ${binDir}adjust/sc_adjustment.py \
         ${(params.containsKey('adjustmentMethod')) ? '--method ' + params.adjustmentMethod : ''} \
         ${(params.containsKey('normalizationVariablesToRegressOut')) ? normalizationVariablesToRegressOutAsArguments : ''} \
         $f \
         "${id}.SC__SCANPY__ADJUSTMENT.${params.off}" 
    """
}
