nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  binDir = ""
}

process SC__SCATER__COMPUTE_QC_STATS {

}

process SC__SCATER__CELL__FILTER {

    container params.sc.scater.container
    
    input:
      tuple val(id), file(f)
    output:
      tuple val(id), file("${id}.SC__SCATER__CELL__FILTER.${params.off}")
    script:
      """
      ${binDir}filter/sc_cell_gene_filtering.py \
        --rdsFile $f \
        --output ${id}.SC__SCATER__CELL_FILTER.${params.off} \
        ${(params.containsKey('nmads')) ? '--nmads ' + params.nmads : ''}
      """

}