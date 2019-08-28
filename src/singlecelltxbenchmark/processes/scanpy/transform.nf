nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

process SC__SCANPY__NORMALIZATION {
  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__NORMALIZATION.${params.off}"
  script:
    """
    python ../../../src/singlecelltxbenchmark/scripts/scanpy/transform/sc_normalization.py \
        ${(params.containsKey('normalizationMethod')) ? '--method ' + params.normalizationMethod : ''} \
        ${(params.containsKey('countsPerCellAfter')) ? '--counts-per-cell-after ' + params.countsPerCellAfter : ''} \
        $f \
        "${getBaseName(f)}.SC__SCANPY__NORMALIZATION.${params.off}"
    """
}

process SC__SCANPY__DATA_TRANSFORMATION {
  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__DATA_TRANSFORMATION.${params.off}"
  script:
    """
    python ../../../src/singlecelltxbenchmark/scripts/scanpy/transform/sc_data_transformation.py \
        ${(params.containsKey('dataTransformationMethod')) ? '--method ' + params.dataTransformationMethod : ''} \
        $f \
        "${getBaseName(f)}.SC__SCANPY__DATA_TRANSFORMATION.${params.off}"
    """
}

process SC__SCANPY__FEATURE_SCALING {
  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__FEATURE_SCALING.${params.off}"
  script:
    """
    python ../../../src/singlecelltxbenchmark/src/scripts/scanpy/transform/sc_feature_scaling.py \
        ${(params.containsKey('featureScalingMthod')) ? '--method ' + params.featureScalingMthod : ''} \
        ${(params.containsKey('featureScalingMaxSD')) ? '--max-sd ' + params.featureScalingMaxSD : ''} \
       $f \
       "${getBaseName(f)}.SC__SCANPY__FEATURE_SCALING.${params.off}"
    """
}
