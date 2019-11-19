nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scater/bin/"
} else {
  binDir = ""
}

process SC__SCATER__CELL_FILTER {

    publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    container params.sc.scater.container
    
    input:
      tuple val(id), file(f)
    output:
      tuple val(id), file("${id}.SC__SCATER__CELL_FILTER.${params.off}")
    script:
      """
      Rscript ${binDir}filter/sc_cell_gene_filtering.R \
        --rdsFile $f \
        --output ${id}.SC__SCATER__CELL_FILTER.${params.off} \
        ${(params.containsKey('nmads')) ? '--nmads ' + params.nmads : ''}
      """

}