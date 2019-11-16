nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  binDir = ""
}

process SC__SCANPY__CLUSTERING {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true
  
  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.SC__SCANPY__CLUSTERING.${params.off}")
  script:
    """
    ${binDir}cluster/sc_clustering.py \
         ${(params.containsKey('clusteringMethod')) ? '--method ' + params.clusteringMethod : ''} \
         ${(params.containsKey('resolution')) ? '--resolution ' + params.resolution : ''} \
         $f \
         "${id}.SC__SCANPY__CLUSTERING.${params.off}"
    """
}
