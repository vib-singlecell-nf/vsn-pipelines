nextflow.preview.dsl=2

include getBaseName from './files.nf'

process SC__H5AD_TO_LOOM {

  container "/ddn1/vol1/staging/leuven/res_00001/software/Scanpy/1.4.3/Scanpy.sif"
  publishDir "${params.outdir}/loom", mode: 'symlink'

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

process SC__H5AD_TO_BARE_LOOM {

  container "/ddn1/vol1/staging/leuven/res_00001/software/Scanpy/1.4.3/Scanpy.sif"
  publishDir "${params.outdir}/loom", mode: 'symlink'

  input:
    file(f)
  output:
    file "${getBaseName(f)}.bare.loom"
  script:
    """
    ${workflow.projectDir}/src/utils/bin/h5ad_to_bare_loom.py \
         $f \
         "${getBaseName(f)}.bare.loom"
    """
}

