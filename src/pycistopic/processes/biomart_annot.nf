nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/pycistopic/bin/" : ""

toolParams = params.tools.pycistopic
processParams = params.tools.pycistopic.biomart_annot

process PYCISTOPIC__BIOMART_ANNOT {

    publishDir "${params.global.outdir}/intermediate/pycistopic/biomart/", mode: 'symlink'
    container toolParams.container
    label 'compute_resources__default'

    output:
        path("biomart_annot.pickle")

    script:
        """
        ${binDir}biomart_annot.py \
            --biomart_dataset_name ${processParams.biomart_dataset_name} \
            --biomart_host ${processParams.biomart_host}
        """
}

