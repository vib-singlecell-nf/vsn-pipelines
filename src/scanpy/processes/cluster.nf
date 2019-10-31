nextflow.preview.dsl=2

include getBaseName from '../../utils/processes/files.nf'

process SC__SCANPY__CLUSTERING {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true
  
  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__CLUSTERING.${params.off}"
  script:
    """
    ${workflow.projectDir}/src/scanpy/bin/cluster/sc_clustering.py \
         ${(params.containsKey('clusteringMethod')) ? '--method ' + params.clusteringMethod : ''} \
         ${(params.containsKey('resolution')) ? '--resolution ' + params.resolution : ''} \
         $f \
         "${getBaseName(f)}.SC__SCANPY__CLUSTERING.${params.off}"
    """
}

