nextflow.enable.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.sinto

process SC__SINTO__FRAGMENTS {

    container toolParams.container
    label 'compute_resources__cpu','compute_resources__24hqueue'

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
            ${processParams.containsKey('temp_dir') && processParams.temp_dir ? '--temp_dir ' + processParams.temp_dir: ''} \
            -p ${task.cpus} \
            -f ${sampleId}.fragments.bed
        """
}


process SC__SINTO__SORT_FRAGMENTS {

    container toolParams.container
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

