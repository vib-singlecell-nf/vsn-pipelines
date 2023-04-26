nextflow.enable.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.sinto

process SINTO__FRAGMENTS {

    container toolParams.container
    label 'compute_resources__sinto__fragments'

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
            ${processParams.containsKey('use_chrom') && processParams.use_chrom ? '--use_chrom ' + processParams.use_chrom: ''} \
            ${processParams.containsKey('min_distance') && processParams.min_distance ? '--min_distance ' + processParams.min_distance: ''} \
            ${processParams.containsKey('max_distance') && processParams.max_distance ? '--max_distance ' + processParams.max_distance: ''} \
            ${processParams.containsKey('chunksize') && processParams.chunksize ? '--chunksize ' + processParams.chunksize: ''} \
            -p ${task.cpus} \
            -f ${sampleId}.fragments.bed
        """
}


process SINTO__SORT_FRAGMENTS {

    container toolParams.container
    label 'compute_resources__sinto__sort_fragments'

    input:
        tuple val(sampleId),
              path(fragments_bed)

    output:
        tuple val(sampleId),
              path("${sampleId}.sinto.fragments.tsv.gz")

    script:
        """
        LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n \
            ${fragments_bed} \
            | bgzip -c \
            > ${sampleId}.sinto.fragments.tsv.gz
        """
}

