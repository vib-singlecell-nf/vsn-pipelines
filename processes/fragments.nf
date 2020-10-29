nextflow.preview.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.sc.atac.sinto

process SC__SINTO__FRAGMENTS {

    container toolParams.container
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId),
              path(bam),
              path(bai)

    output:
        tuple val(sampleId),
              path("${sampleId}.fragments.bed")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.fragments)
        processParams = sampleParams.local
        """
        sinto fragments \
            -b ${bam} \
            -m ${processParams.min_mapq} \
            ${processParams.containsKey('barcodetag') && processParams.barcodetag ? '--barcodetag ' + processParams.barcodetag: ''} \
            ${processParams.containsKey('barcode_regex') && processParams.barcode_regex ? '--barcode_regex ' + processParams.barcode_regex: ''} \
            ${processParams.containsKey('min_distance') && processParams.min_distance ? '--min_distance ' + processParams.min_distance: ''} \
            ${processParams.containsKey('max_distance') && processParams.max_distance ? '--max_distance ' + processParams.max_distance: ''} \
            ${processParams.containsKey('chunksize') && processParams.chunksize ? '--chunksize ' + processParams.chunksize: ''} \
            -p ${task.cpus} \
            -f ${sampleId}.fragments.bed
        """
}


process SC__SINTO__SORT_FRAGMENTS {

    container toolParams.container
    publishDir "${params.global.outdir}/fragments/", mode: 'link'
    label 'compute_resources__mem'

    input:
        tuple val(sampleId),
              path(fragments_bed)

    output:
        tuple val(sampleId),
              path("${sampleId}.sinto.fragments.tsv.gz")

    script:
        """
        sort -k1,1 -k2,2n \
            ${fragments_bed} \
            | bgzip -c \
            > ${sampleId}.sinto.fragments.tsv.gz
        """
}

