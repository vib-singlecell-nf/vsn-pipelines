nextflow.preview.dsl=2

include getBaseName from '../../utils/processes/files.nf'

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  binDir = ""
}

process SC__SCANPY__COMPUTE_QC_STATS {

  container params.sc.scanpy.container

  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__COMPUTE_QC_STATS.${params.off}"
  script:
    """
    ${binDir}filter/sc_cell_gene_filtering.py \
      compute \
      $f \
      ${getBaseName(f)}.SC__SCANPY__COMPUTE_QC_STATS.${params.off} \
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
    publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true

    input:
        file(f)
    output:
        file "${getBaseName(f)}.SC__SCANPY__GENE_FILTER.${params.off}"
    script:
    """
    ${binDir}filter/sc_cell_gene_filtering.py \
        genefilter \
        $f \
        ${getBaseName(f)}.SC__SCANPY__GENE_FILTER.${params.off} \
        ${(params.containsKey('geneFilterMinNCells')) ? '--min-number-cells ' + params.geneFilterMinNCells : ''} \
    """
}


process SC__SCANPY__CELL_FILTER {

    container params.sc.scanpy.container
    publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true

    input:
        file(f)
    output:
        file "${getBaseName(f)}.SC__SCANPY__CELL_FILTER.${params.off}"
    script:
    """
    ${binDir}filter/sc_cell_gene_filtering.py \
        cellfilter \
        $f \
        ${getBaseName(f)}.SC__SCANPY__CELL_FILTER.${params.off} \
        ${(params.containsKey('cellFilterMinNCounts')) ? '--min-n-counts ' + params.cellFilterMinNCounts : ''} \
        ${(params.containsKey('cellFilterMaxNCounts')) ? '--max-n-counts ' + params.cellFilterMaxNCounts : ''} \
        ${(params.containsKey('cellFilterMinNGenes')) ? '--min-n-genes ' + params.cellFilterMinNGenes : ''} \
        ${(params.containsKey('cellFilterMaxNGenes')) ? '--max-n-genes ' + params.cellFilterMaxNGenes : ''} \
        ${(params.containsKey('cellFilterMaxPercentMito')) ? '--max-percent-mito ' + params.cellFilterMaxPercentMito : ''}
    """
}
