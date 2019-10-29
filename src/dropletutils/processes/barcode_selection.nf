nextflow.preview.dsl=2

process SC__DROPLET_UTILS__BARCODE_SELECTION {
    
    container params.sc.dropletutils.container
    publishDir "03.count", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(readCounts)
    output:
		tuple val(sample), val("knee"), file("*.selected_cell_barcodes_by_knee.txt"), emit: selectedCellBarcodesByKnee
        tuple val(sample), val("inflection"), file("*.selected_cell_barcodes_by_inflection.txt"), emit: selectedCellBarcodesByInflection
        tuple file("*.barcode_rank_vs_total_umi_plot.png"), emit: plot
    script:
        """
        Rscript ${workflow.projectDir}/src/dropletutils/bin/barcode_selection.R \
            ${readCounts} \
            ${sample}
        """    
}