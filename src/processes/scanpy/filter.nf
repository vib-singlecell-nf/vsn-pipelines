nextflow.preview.dsl=2

include getBaseName from '../../utils/files.nf'

params.geneFilterMinNCells = -1

process SC__SCANPY__GENE_FILTER {
  input:
    file(f)
  output:
    file "${getBaseName(f).get()}.SC__SCANPY__GENE_FILTER.${params.off}"
  script:
    """
    python ../../../src/scripts/scanpy/filter/sc_gene_filter.py \
        --min-number-cells $params.geneFilterMinNCells \
       $f \
       "${getBaseName(f).get()}.SC__SCANPY__GENE_FILTER.${params.off}"
    """
}

params.cellFilterMinNCounts = -1
params.cellFilterMaxNCounts = -1
params.cellFilterMinNGenes = -1
params.cellFilterMaxNGenes = -1
params.cellFilterMaxPercentMito = -1

process SC__SCANPY__CELL_FILTER {
  input:
    file(f)
  output:
    file "${getBaseName(f).get()}.SC__SCANPY__CELL_FILTER.${params.off}"
  script:
    """
    python ../../../src/scripts/scanpy/filter/sc_cell_filter.py \
        --min-n-counts $params.cellFilterMinNCounts \
        --max-n-counts $params.cellFilterMaxNCounts \
        --min-n-genes $params.cellFilterMinNGenes \
        --max-n-genes $params.cellFilterMaxNGenes \
        --max-percent-mito $params.cellFilterMaxPercentMito \
        $f \
        "${getBaseName(f).get()}.SC__SCANPY__CELL_FILTER.${params.off}"
    """
}

process SC__SCANPY__FILTER_QC_REPORT {
  input:
    file(f)
  output:
    file "${getBaseName(f).get()}.SC__SCANPY__FILTER_QC_REPORT.ipynb"
  script:
    """
    papermill ../../../src/scripts/scanpy/filter/sc_filter_qc_report.ipynb \
        ${getBaseName(f).get()}.SC__SCANPY__FILTER_QC_REPORT.ipynb \
        -p FILE $f
    """
}
