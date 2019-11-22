nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  binDir = ""
}

process SC__SCANPY__COMPUTE_QC_STATS {

  container params.sc.scanpy.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"

  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.SC__SCANPY__COMPUTE_QC_STATS.${params.off}")
  script:
    """
    ${binDir}filter/sc_cell_gene_filtering.py \
      compute \
      $f \
      ${id}.SC__SCANPY__COMPUTE_QC_STATS.${params.off} \
      ${(params.containsKey('cellFilterMinNCounts')) ? '--min-n-counts ' + params.cellFilterMinNCounts : ''} \
      ${(params.containsKey('cellFilterMaxNCounts')) ? '--max-n-counts ' + params.cellFilterMaxNCounts : ''} \
      ${(params.containsKey('cellFilterMinNGenes')) ? '--min-n-genes ' + params.cellFilterMinNGenes : ''} \
      ${(params.containsKey('cellFilterMaxNGenes')) ? '--max-n-genes ' + params.cellFilterMaxNGenes : ''} \
      ${(params.containsKey('cellFilterMaxPercentMito')) ? '--max-percent-mito ' + params.cellFilterMaxPercentMito : ''} \
      ${(params.containsKey('geneFilterMinNCells')) ? '--min-number-cells ' + params.geneFilterMinNCells : ''}
    """
}


process SC__SCANPY__GENE_FILTER {

    container params.sc.scanpy.container
    clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true

    input:
        tuple val(id), file(f)
    output:
        tuple val(id), file("${id}.SC__SCANPY__GENE_FILTER.${params.off}")
    script:
    """
    ${binDir}filter/sc_cell_gene_filtering.py \
        genefilter \
        $f \
        ${id}.SC__SCANPY__GENE_FILTER.${params.off} \
        ${(params.containsKey('geneFilterMinNCells')) ? '--min-number-cells ' + params.geneFilterMinNCells : ''}
    """
}


process SC__SCANPY__CELL_FILTER {

    container params.sc.scanpy.container
    clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true

    input:
        tuple val(id), file(f)
    output:
        tuple val(id), file("${id}.SC__SCANPY__CELL_FILTER.${params.off}")
    script:
    """
    ${binDir}filter/sc_cell_gene_filtering.py \
        cellfilter \
        $f \
        ${id}.SC__SCANPY__CELL_FILTER.${params.off} \
        ${(params.containsKey('cellFilterMinNCounts')) ? '--min-n-counts ' + params.cellFilterMinNCounts : ''} \
        ${(params.containsKey('cellFilterMaxNCounts')) ? '--max-n-counts ' + params.cellFilterMaxNCounts : ''} \
        ${(params.containsKey('cellFilterMinNGenes')) ? '--min-n-genes ' + params.cellFilterMinNGenes : ''} \
        ${(params.containsKey('cellFilterMaxNGenes')) ? '--max-n-genes ' + params.cellFilterMaxNGenes : ''} \
        ${(params.containsKey('cellFilterMaxPercentMito')) ? '--max-percent-mito ' + params.cellFilterMaxPercentMito : ''}
    """
}
