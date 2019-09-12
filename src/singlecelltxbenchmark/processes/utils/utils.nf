nextflow.preview.dsl=2

// Specifiy the input and out format files
// 'if' not accepted as param
// 'in' not allwed as variable
params.project_name = ''
params.iff = ''
params.off = ''

process SC__FILE_CONVERTER {

  publishDir "${params.outdir}/data", mode: 'symlink'

  input:
    set id, file(f)
  output:
    file "${id}.SC__FILE_CONVERTER.${params.off}"
  script:
    """
    python $params.baseFilePath/src/singlecelltxbenchmark/scripts/utils/sc_file_converter.py \
       --input-format $params.iff \
       --output-format $params.off $f "${id}.SC__FILE_CONVERTER.${params.off}"
    """
}

process SC__FILE_CONVERTER_HELP {
  output:
    stdout()
  script:
    """
    python $params.baseFilePath/src/singlecelltxbenchmark/scripts/utils/sc_file_converter.py -h | awk '/-h/{y=1;next}y'
    """
}

process SC__FILE_CONCATENATOR() {

  publishDir "${params.outdir}/data", mode: 'symlink'

  input:
    file(f)
  output:
    file "${params.project_name}.SC__FILE_CONCATENATOR.${params.off}"
  script:
    """
    $params.baseFilePath/src/singlecelltxbenchmark/scripts/utils/sc_file_concatenator.py \
      --file-format $params.off \
      ${(params.containsKey('join')) ? '--join ' + params.join : ''} \
      --output "${params.project_name}.SC__FILE_CONCATENATOR.${params.off}" $f
    """
}

include getBaseName from '../../utils/files.nf'

process SC__FILE_ANNOTATOR() {

  publishDir "${params.outdir}/data", mode: 'symlink'

  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__FILE_ANNOTATOR.${params.off}"
  script:
    """
    $params.baseFilePath/src/singlecelltxbenchmark/scripts/utils/sc_file_annotator.py \
      ${(params.containsKey('type')) ? '--type ' + params.type : ''} \
      ${(params.containsKey('metaDataFilePath')) ? '--meta-data-file-path ' + params.metaDataFilePath : ''} \
      $f \
      "${getBaseName(f)}.SC__FILE_ANNOTATOR.${params.off}"
    """
}