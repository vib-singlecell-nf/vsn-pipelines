nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

params.normalizationMethod = 'cpx'
params.normalizationNumberReadsPerCellFactor = -1

process SC__SCANPY__NORMALIZATION {
  input:
    file(f)
  output:
    file "${getBaseName(f).get()}.SC__SCANPY__NORMALIZATION.${params.off}"
  script:
    """
    python ../../../src/scripts/scanpy/transform/sc_normalization.py \
        --method $params.normalizationMethod \
        --number-reads-per-cell-factor $params.normalizationNumberReadsPerCellFactor \
        $f \
        "${getBaseName(f).get()}.SC__SCANPY__NORMALIZATION.${params.off}"
    """
}

params.dataTransformationMethod = 'log1p'

process SC__SCANPY__DATA_TRANSFORMATION {
  input:
    file(f)
  output:
    file "${getBaseName(f).get()}.SC__SCANPY__DATA_TRANSFORMATION.${params.off}"
  script:
    """
    python ../../../src/scripts/scanpy/transform/sc_data_transformation.py \
        --method $params.dataTransformationMethod \
        $f \
        "${getBaseName(f).get()}.SC__SCANPY__DATA_TRANSFORMATION.${params.off}"
    """
}

params.featureScalingMthod = 'zscore_scale'
params.featureScalingMaxSD = -1

process SC__SCANPY__FEATURE_SCALING {
  input:
    file(f)
  output:
    file "${getBaseName(f).get()}.SC__SCANPY__FEATURE_SCALING.${params.off}"
  script:
    """
    python ../../../src/scripts/scanpy/transform/sc_feature_scaling.py \
        --method $params.featureScalingMthod \
        --max-sd $params.featureScalingMaxSD \
       $f \
       "${getBaseName(f).get()}.SC__SCANPY__FEATURE_SCALING.${params.off}"
    """
}
