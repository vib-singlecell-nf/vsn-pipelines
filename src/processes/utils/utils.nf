nextflow.preview.dsl=2

// Specifiy the input and out format files
// 'if' not accepted as param
// 'in' not allwed as variable
params.project_name = ''
params.iff = ''
params.off = ''

process SC__FILE_CONVERTER {
  input:
    set id, file(f)
  output:
    file "${id}.SC__FILE_CONVERTER.${params.off}"
  script:
    """
    python ../../../src/scripts/utils/sc_file_converter.py --input-format $params.iff --output-format $params.off $f "${id}.SC__FILE_CONVERTER.${params.off}"
    """
}

process SC__FILE_CONCATENATOR() {
  input:
    file(f)
  output:
    file "${params.project_name}.SC__FILE_CONCATENATOR.${params.off}"
  script:
    """
    python ../../../src/scripts/utils/sc_file_concatenator.py \
      --file-format $params.off \
      ${(params.containsKey('join')) ? '--join ' + params.join : ''} \
      --output "${params.project_name}.SC__FILE_CONCATENATOR.${params.off}" $f
    """
}