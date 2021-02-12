nextflow.enable.dsl=2

process SC__DROPLET_UTILS__BARCODE_SELECTION {
    
    container params.getToolParams("dropletutils")..container
    publishDir "03.count", mode: 'symlink'
    label 'compute_resources__default'

    input:
        tuple val(sample), path(readCounts)
    
    output:
        tuple val(sample), val("knee"), path("*.selected_cell_barcodes_by_knee.txt"), emit: selectedCellBarcodesByKnee
        tuple val(sample), val("inflection"), path("*.selected_cell_barcodes_by_inflection.txt"), emit: selectedCellBarcodesByInflection
        tuple file("*.barcode_rank_vs_total_umi_plot.png"), emit: plot
    
    script:
        """
        Rscript ${workflow.projectDir}/src/dropletutils/bin/barcode_selection.R \
            ${readCounts} \
            ${sample}
        """

}
