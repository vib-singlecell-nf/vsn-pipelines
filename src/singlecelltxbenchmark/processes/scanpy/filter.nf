nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

process SC__SCANPY__GENE_FILTER {

  publishDir "${params.outdir}/data", mode: 'symlink'

  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__GENE_FILTER.${params.off}"
  script:
    """
    sc_gene_filter.py \
      ${(params.containsKey('geneFilterMinNCells')) ? '--min-number-cells ' + params.geneFilterMinNCells : ''} \
      $f \
      "${getBaseName(f)}.SC__SCANPY__GENE_FILTER.${params.off}"
    """
}

process SC__SCANPY__CELL_FILTER {

  publishDir "${params.outdir}/data", mode: 'symlink'

  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__CELL_FILTER.${params.off}"
  script:
    """
    sc_cell_filter.py \
        ${(params.containsKey('cellFilterMinNCounts')) ? '--min-n-counts ' + params.cellFilterMinNCounts : ''} \
        ${(params.containsKey('cellFilterMaxNCounts')) ? '--max-n-counts ' + params.cellFilterMaxNCounts : ''} \
        ${(params.containsKey('cellFilterMinNGenes')) ? '--min-n-genes ' + params.cellFilterMinNGenes : ''} \
        ${(params.containsKey('cellFilterMaxNGenes')) ? '--max-n-genes ' + params.cellFilterMaxNGenes : ''} \
        ${(params.containsKey('cellFilterMaxPercentMito')) ? '--max-percent-mito ' + params.cellFilterMaxPercentMito : ''} \
        $f \
        "${getBaseName(f)}.SC__SCANPY__CELL_FILTER.${params.off}"
    """
}

process SC__SCANPY__FILTER_QC_REPORT {

  publishDir "${params.outdir}/notebooks", mode: 'symlink'

  input:
    file(f)
  output:
    file "${getBaseName(f)}.SC__SCANPY__FILTER_QC_REPORT.ipynb"
  script:
    """
    papermill $params.baseFilePath/src/singlecelltxbenchmark/scripts/scanpy/filter/sc_filter_qc_report.ipynb \
        ${getBaseName(f)}.SC__SCANPY__FILTER_QC_REPORT.ipynb \
        -p FILE $f
    """
}
