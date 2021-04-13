nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.bap

process SC__BAP__BIORAD_DEBARCODE {

    container toolParams.container
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sampleId),
              path(fastq_PE1),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}*1.fastq.gz"),
              path("${sampleId}*2.fastq.gz"),
              path("${sampleId}-parse.sumstats.log")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.biorad_debarcode)
        processParams = sampleParams.local
        """
        bap-barcode ${processParams.protocol} \
            -a ${fastq_PE1} \
            -b ${fastq_PE2} \
            --ncores ${task.cpus} \
            --nmismatches ${processParams.nmismatches} \
            --output ${sampleId} 
        """
}


process SC__BAP__MERGE_FASTQS {

    container toolParams.container
    label 'compute_resources__default'

    input:
        tuple val(sampleId),
              path(fastq_PE1),
              path(fastq_PE2),
              path(log)

    output:
        tuple val(sampleId),
              path("${sampleId}_dex_R1.fastq.gz"),
              path("${sampleId}_dex_R2.fastq.gz")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        """
        cat ${fastq_PE1} > ${sampleId}_dex_R1.fastq.gz
        cat ${fastq_PE2} > ${sampleId}_dex_R2.fastq.gz
        """
}

