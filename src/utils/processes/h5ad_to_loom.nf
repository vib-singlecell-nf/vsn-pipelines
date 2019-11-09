nextflow.preview.dsl=2

include getBaseName from './files.nf'

process SC__H5AD_TO_LOOM {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/loom", mode: 'link', overwrite: true

  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__H5AD_TO_LOOM.loom" 
  script:
    """
    ${workflow.projectDir}/src/utils/bin/h5ad_to_loom.py \
         $f \
         "${getBaseName(f)}.SC__H5AD_TO_LOOM.loom" 
    """
}

process SC__H5AD_TO_FILTERED_LOOM {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/loom", mode: 'link', overwrite: true

  input:
    file(f)
  output:
    file "${getBaseName(f)}.filtered.loom"
  script:
    """
    ${workflow.projectDir}/src/utils/bin/h5ad_to_filtered_loom.py \
         $f \
         "${getBaseName(f)}.filtered.loom"
    """
}

