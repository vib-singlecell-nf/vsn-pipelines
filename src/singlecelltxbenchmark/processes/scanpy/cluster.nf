nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

process SC__SCANPY__CLUSTERING {

  publishDir "${params.outdir}/data", mode: 'symlink'
  
  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__CLUSTERING.${params.off}"
  script:
    """
    sc_clustering.py \
         ${(params.containsKey('clusteringMethod')) ? '--method ' + params.clusteringMethod : ''} \
         ${(params.containsKey('resolution')) ? '--resolution ' + params.resolution : ''} \
         $f \
         "${getBaseName(f)}.SC__SCANPY__CLUSTERING.${params.off}"
    """
}
