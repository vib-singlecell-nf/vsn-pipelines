nextflow.preview.dsl=2

include getBaseName from '../../utils/processes/files.nf'

process SC__SCANPY__MARKER_GENES {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true
  
  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__MARKER_GENES.${params.off}"
  script:
    """
    ${workflow.projectDir}/src/scanpy/bin/cluster/sc_marker_genes.py \
         ${(params.containsKey('method')) ? '--method ' + params.method : ''} \
         ${(params.containsKey('groupby')) ? '--groupby ' + params.groupby : ''} \
         ${(params.containsKey('ngenes')) ? '--ngenes ' + params.ngenes : ''} \
         $f \
         "${getBaseName(f)}.SC__SCANPY__MARKER_GENES.${params.off}"
    """
}

