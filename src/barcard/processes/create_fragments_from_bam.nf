nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""


process CREATE_FRAGMENTS_FROM_BAM {
    container params.tools.barcard.container
    label 'compute_resources__barcard__create_fragments_from_bam'
    publishDir "${params.global.outdir}/data/fragments", mode: 'copy'

    input:
        tuple val(sampleId),
              path(bam),
              path(bai)

    output:
        tuple val(sampleId),
              path("${sampleId}.fragments.raw.tsv.gz"),
              path("${sampleId}.fragments.raw.tsv.gz.tbi")

    script:
        //def sampleParams = params.parseConfig(sampleId, params.global)
        //processParams = sampleParams.local
        """
        set -euo pipefail
        mkdir ./tmp
        export TMPDIR=./tmp
        create_fragments_file \
            "${bam}" \
          | coreutils sort --parallel=8 -S 16G -k 1,1V -k 2,2n -k 3,3n -k 4,4 \
          | uniq -c \
          | mawk -v 'OFS=\t' '{ print \$2, \$3, \$4, \$5, \$1 }' \
          | bgzip -@ 4 -c /dev/stdin \
          > ${sampleId}.fragments.raw.tsv.gz

        tabix -p bed ${sampleId}.fragments.raw.tsv.gz
        """
}
