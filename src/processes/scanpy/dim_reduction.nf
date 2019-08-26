nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

params.dimReductionMethod = 'PCA'

process SC__SCANPY__DIM_REDUCTION {
  input:
    file(f)
  output:
    file "${getBaseName(f).get()}.SC__SCANPY__DIM_REDUCTION.${params.off}"
  script:
    """
    python ../../../src/scripts/scanpy/dim_reduction/sc_dim_reduction.py \
         --method $params.dimReductionMethod \
         ${(params.containsKey('svdSolver')) ? '--svd-solver ' + params.svdSolver : ''} \
         ${(params.containsKey('nNeighbors')) ? '--n-neighbors ' + params.nNeighbors : ''} \
         ${(params.containsKey('nComps')) ? '--n-comps ' + params.nComps : ''} \
         ${(params.containsKey('nPcs')) ? '--n-pcs ' + params.nPcs : ''} \
         ${(params.containsKey('nJobs')) ? '--n-jobs ' + params.nJobs : ''} \
         $f \
         "${getBaseName(f).get()}.SC__SCANPY__DIM_REDUCTION.${params.off}"
    """
}
