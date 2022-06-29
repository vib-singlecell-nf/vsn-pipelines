nextflow.enable.dsl=2

process SC__SEURAT__GENERATE_REPORT {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/markdowns", mode: 'link', overwrite: true
    label 'compute_resources__report'

    input:
        file(rmd)
        tuple val(sampleId), path(seurat_rds)
        val(reportTitle)

    output:
        tuple val(sampleId), path("${sampleId}.${reportTitle}.Rmd")

    script:
        """
        cat ${rmd} | sed 's/FILE/${seurat_rds}/' > ${sampleId}.${reportTitle}.Rmd
        """
}

process SC__SEURAT__GENERATE_DUAL_INPUT_REPORT {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/markdowns", mode: 'link', overwrite: true
    label 'compute_resources__report'

    input:
        file(rmd)
        tuple val(sampleId), path(seurat_rds1)
        tuple val(sampleId2), path(seurat_rds2)
        val(reportTitle)

    output:
        tuple val(sampleId), path("${sampleId}.${reportTitle}.Rmd")

    script:
        """
        cat ${rmd} | sed 's/FILE1/${seurat_rds1}/' | sed 's/FILE2/${seurat_rds2}/' > ${sampleId}.${reportTitle}.Rmd
        """
}

process SC__SEURAT__REPORT_TO_HTML {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/reports", mode: 'link', overwrite: true
    label 'compute_resources__report'

    input:
        tuple val(sampleId), path(rmd)
    
    output:
        file("*.html")
    
    // Running rmarkdown::render on a symlinked file, seems to create the rendered html in the original location.
    // To get the output html in the correct folder, we need to set the output_dir parameter to the current working directory.
    shell:
        '''
        R -e "rmarkdown::render(
            input = '!{rmd}',
            output_dir = '$(pwd)'
        )"
        '''
}
