nextflow.enable.dsl=2

process SC__SEURAT__GENERATE_REPORT {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/markdowns/intermediate", mode: 'link', overwrite: true
    label 'compute_resources__report'

    input:
        file(rmd)
        tuple val(sampleId), path(seurat_rds)
        val(reportTitle)

    output:
        tuple val(sampleId), path("${sampleId}.${reportTitle}.rmd")

    script:
        """
        cat ${rmd} | sed 's/FILE/${seurat_rds}/' > ${sampleId}.${reportTitle}.rmd
        """
}

process SC__SEURAT__REPORT_TO_HTML {

    // container params.tools.seurat.container
    publishDir "${params.global.outdir}/markdowns/intermediate", mode: 'link', overwrite: true
    label 'compute_resources__report'

    input:
        tuple val(sampleId), path(rmd)
        tuple val(samleId), path(rds)
        val(reportTitle)
    
    output:
        file("${sampleId}.${reportTitle}.html")
    
    // Running rmarkdown::render on a symlinked file, seems to create the rendered html in the original location.
    // To get the output html in the correct folder, we need to set the output_dir parameter to the current working directory.
    shell:
        '''
        R -e "rmarkdown::render(
            input = '!{rmd}',
            output_format = 'html_document',
            output_file = '!{sampleId}.!{reportTitle}.html',
            output_dir = '$(pwd)'
        )"
        '''
}
