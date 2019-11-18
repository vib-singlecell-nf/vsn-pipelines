nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  binDir = ""
}

process SC__SCANPY__MARKER_GENES {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true
  
  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.SC__SCANPY__MARKER_GENES.${params.off}")
  script:
    """
    ${binDir}cluster/sc_marker_genes.py \
         ${(params.containsKey('method')) ? '--method ' + params.method : ''} \
         ${(params.containsKey('groupby')) ? '--groupby ' + params.groupby : ''} \
         ${(params.containsKey('ngenes')) ? '--ngenes ' + params.ngenes : ''} \
         $f \
         "${id}.SC__SCANPY__MARKER_GENES.${params.off}"
    """
}

