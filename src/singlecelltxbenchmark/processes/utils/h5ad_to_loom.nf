nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

process SC__H5AD_TO_LOOM {

  publishDir "${params.outdir}/loom", mode: 'symlink'
  container "/ddn1/vol1/staging/leuven/res_00001/software/Scanpy/1.4.3/Scanpy.sif"

  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__H5AD_TO_LOOM.loom" 
  script:
    """
    h5ad_to_loom.py \
         $f \
         "${getBaseName(f)}.SC__H5AD_TO_LOOM.loom" 
    """
}
