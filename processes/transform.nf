nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  binDir = ""
}

process SC__SCANPY__NORMALIZATION {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.SC__SCANPY__NORMALIZATION.${params.off}")
  script:
    """
    ${binDir}transform/sc_normalization.py \
        ${(params.containsKey('normalizationMethod')) ? '--method ' + params.normalizationMethod : ''} \
        ${(params.containsKey('countsPerCellAfter')) ? '--counts-per-cell-after ' + params.countsPerCellAfter : ''} \
        $f \
        "${id}.SC__SCANPY__NORMALIZATION.${params.off}"
    """
}

process SC__SCANPY__DATA_TRANSFORMATION {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.SC__SCANPY__DATA_TRANSFORMATION.${params.off}")
  script:
    """
    ${binDir}transform/sc_data_transformation.py \
        ${(params.containsKey('dataTransformationMethod')) ? '--method ' + params.dataTransformationMethod : ''} \
        $f \
        "${id}.SC__SCANPY__DATA_TRANSFORMATION.${params.off}"
    """
}

process SC__SCANPY__FEATURE_SCALING {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.SC__SCANPY__FEATURE_SCALING.${params.off}")
  script:
    """
    ${binDir}transform/sc_feature_scaling.py \
        ${(params.containsKey('featureScalingMthod')) ? '--method ' + params.featureScalingMthod : ''} \
        ${(params.containsKey('featureScalingMaxSD')) ? '--max-sd ' + params.featureScalingMaxSD : ''} \
       $f \
       "${id}.SC__SCANPY__FEATURE_SCALING.${params.off}"
    """
}
