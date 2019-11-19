nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  binDir = ""
}

process SC__SCANPY__DIM_REDUCTION {

  container params.sc.scanpy.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.SC__SCANPY__DIM_REDUCTION_${method}.${params.off}")
  script:
    method = params.dimReductionMethod.replaceAll('-','').toUpperCase()
    """
    ${binDir}dim_reduction/sc_dim_reduction.py \
         ${(params.containsKey('dimReductionMethod')) ? '--method ' + params.dimReductionMethod : ''} \
         ${(params.containsKey('svdSolver')) ? '--svd-solver ' + params.svdSolver : ''} \
         ${(params.containsKey('nNeighbors')) ? '--n-neighbors ' + params.nNeighbors : ''} \
         ${(params.containsKey('nComps')) ? '--n-comps ' + params.nComps : ''} \
         ${(params.containsKey('nPcs')) ? '--n-pcs ' + params.nPcs : ''} \
         ${(params.containsKey('nJobs')) ? '--n-jobs ' + params.nJobs : ''} \
         ${(params.containsKey('useFastTsne') && !params.useFastTsne) ? '' : '--use-fast-tsne'} \
         $f \
         "${id}.SC__SCANPY__DIM_REDUCTION_${method}.${params.off}"
    """
}

