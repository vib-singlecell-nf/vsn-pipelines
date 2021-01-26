nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.bwamaptools

process SC__BWAMAPTOOLS__ADD_BARCODE_TAG {

    container toolParams.container
    label 'compute_resources__default','compute_resources__24hqueue'
    // todo: add storeDir instead of publishDir

    input:
        tuple val(sampleId),
              path(bam)

    output:
        tuple val(sampleId),
              path("${sampleId}.bwa.possorted.bam")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.add_barcode_as_tag)
        processParams = sampleParams.local
        """
        samtools view -h \
            ${bam} \
            | mawk '/^@/ {print;next} {N=split(\$1,n,"${processParams.delimiter_to_split_qname}");print \$0 "\t${processParams.tag}:Z:" n[${processParams.position_of_barcode_in_qname}]}' \
            | samtools view -bS - -o ${sampleId}.bwa.possorted.bam
        """
}

/*
original mawk command:
samtools view -h possorted.bam \
    | mawk '/^@/ {print;next} {N=split($1,n,":");print $0 "\tCR:Z:" n[1]}' \
    | samtools view -bS - -o possorted.bc.bam
*/

