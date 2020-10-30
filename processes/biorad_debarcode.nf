nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.sc.atac.bap

process SC__BAP__BIORAD_DEBARCODE {

    container toolParams.container
    //publishDir "${params.global.outdir}/bap", mode: 'symlink'
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId),
              path(fastq_PE1),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}*1.fastq.gz"),
              path("${sampleId}*2.fastq.gz")
              //path("*/*parse.sumstats.log")

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
    publishDir "${params.global.outdir}/fastq/barcode_demultiplexed", mode: 'symlink'
    label 'compute_resources__default'

    input:
        tuple val(sampleId),
              path(fastq_PE1),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}.R1.dex.fastq.gz"),
              path("${sampleId}.R2.dex.fastq.gz")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
		processParams = sampleParams.local
        """
        cat ${fastq_PE1} > ${sampleId}.R1.dex.fastq.gz
        cat ${fastq_PE2} > ${sampleId}.R2.dex.fastq.gz
        """
}

